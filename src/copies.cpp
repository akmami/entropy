// copies: the repeat structure of a genome from its suffix array and LCP array (libdivsufsort).
//
//   copies genome.fa.gz k radii ks      e.g.  copies genomes/real/ecoli.fa.gz 21 7,10 15,17,19,21,25,31
//
// radii are the context radii of the sketching schemes (src/ue settings prints them), ks the k-mer
// lengths of ue's k-sweep. Writes genome.copies.tsv next to the FASTA, which ue real reads.
//
// The text G is the records joined by 0, every non-ACGT base replaced by 0; matches never cross
// a 0, so nothing spans an N or a record end. One strand only, as in ue: a k-mer and its reverse
// complement are different k-mers. For every k-mer position i of G (k bases, no N):
//
//   c(i)     occurrences of the k-mer in G: ue's m_G
//   m_r(i)   occurrences of the (k + 2r)-mer centred on it: the copies whose whole context of
//            radius r is identical, which a scheme of context radius r selects together
//   R(i)     bases from i on which all c occurrences agree (>= k), from the suffix array of G
//   Lext(i)  bases to the left of i on which all c occurrences agree, from the suffix array of
//            G reversed (not complemented), where the k-mer at i starts at g - i - k
//   extent   Lext + R: the stretch around i that all c copies share
//
// Copies and families are counted through the positions: a family of c identical copies of
// length L has c (L - k + 1) positions, each with count c and extent L, so a position counts as
// 1 / (extent - k + 1) copies and 1 / (c (extent - k + 1)) families. (Counting the left ends of
// the extents instead would also count the chance sub-families at copy edges, where some copies
// happen to share flanking bases: about twice the planted copies on a simulated genome.)
//
// Output: '#key<TAB>value' header lines, then one table row per line, its first column naming
// the table:
//
//   S  c  positions           positions per count c (the position spectrum)
//   P  r  c  m  positions  s_0 .. s_{r-1}
//                             positions with c >= 2 and context count m = m_r, per radius r; s_e:
//                             (position, copy) pairs whose contexts agree on exactly e bases on
//                             either side (and differ at the next base on at least one side)
//   E  extent  positions      extents of the positions with c >= 2
//   F  c  families            families per copy number c
//   C  length  copies         copies per length (the extent)
//   K  k'  c  positions       positions per count for the k-sweep's k-mer lengths k': one pass
//                             over the LCP array per k'
//
// In P, E and C, values above 1024 are rounded to a 1% geometric grid (S and F are exact).
// F and C, and the #copies and #families header lines, are fractional sums.
//
// Memory peaks while an LCP array is built: the text (1 byte per base), the suffix array and Phi
// (8 bytes each), 4 bytes of LCP and 4 bytes of left extensions: about 25 bytes per base (some
// 80 GB for the human genome).
#include <divsufsort64.h>
#include <sys/resource.h>
#include <zlib.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>

typedef saidx64_t idx;  // suffix array index; 64 bits, as genomes over 2^31 bases need

static double t_start;

static double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static double peak_gb() {
    struct rusage u;
    getrusage(RUSAGE_SELF, &u);
#ifdef __APPLE__
    return u.ru_maxrss / 1e9;  // bytes
#else
    return u.ru_maxrss / 1e6;  // kilobytes
#endif
}

static void log(const std::string& what) {
    fprintf(stderr, "[copies %7.1fs %6.1f GB] %s\n", now() - t_start, peak_gb(), what.c_str());
}

// ---------------------------------------------------------------------------
// FASTA.gz -> codes A,C,G,T = 1..4, other bases 0; records joined by 0
// ---------------------------------------------------------------------------
struct Genome {
    std::vector<uint8_t> text;
    long records = 0;
    uint64_t acgt = 0;
};

static Genome read_fasta(const char* path) {
    uint8_t code[256] = {0};
    code['A'] = code['a'] = 1, code['C'] = code['c'] = 2, code['G'] = code['g'] = 3, code['T'] = code['t'] = 4;
    Genome g;
    gzFile f = gzopen(path, "rb");
    if (!f) perror(path), exit(1);
    static char buf[1 << 20];
    bool header = false, line_start = true;
    for (int n; (n = gzread(f, buf, sizeof buf)) > 0;) {
        for (int i = 0; i < n; i++) {
            unsigned char ch = buf[i];
            if (ch == '\n') {
                header = false, line_start = true;
                continue;
            }
            if (line_start && ch == '>') {
                if (!g.text.empty() && g.text.back() != 0) g.text.push_back(0);  // record end
                header = true, g.records++;
            } else if (!header && ch > ' ') {
                g.text.push_back(code[ch]);
                g.acgt += code[ch] != 0;
            }
            line_start = false;
        }
    }
    gzclose(f);
    while (!g.text.empty() && g.text.back() == 0) g.text.pop_back();
    return g;
}

// ---------------------------------------------------------------------------
// suffix array (divsufsort) and LCP array (Kasai via Phi); lcp[r] = common prefix of the
// suffixes SA[r-1] and SA[r], stopping at a 0, so lcp >= x means x shared bases without N
// ---------------------------------------------------------------------------
static void suffix_lcp(const std::vector<uint8_t>& T, std::vector<idx>& SA, std::vector<uint32_t>& lcp) {
    const idx n = (idx)T.size();
    SA.resize(n);
    divsufsort64(T.data(), SA.data(), n);
    log("suffix array");
    std::vector<uint32_t> plcp(n);  // by text position
    {
        std::vector<idx> phi(n);  // phi[SA[r]] = SA[r-1]
        phi[SA[0]] = -1;
        for (idx r = 1; r < n; r++) phi[SA[r]] = SA[r - 1];
        for (idx i = 0, h = 0; i < n; i++) {
            idx j = phi[i];
            if (j < 0) {
                plcp[i] = 0, h = 0;
                continue;
            }
            while (i + h < n && j + h < n && T[i + h] && T[i + h] == T[j + h]) h++;
            plcp[i] = (uint32_t)h;
            if (h > 0) h--;
        }
    }
    lcp.resize(n);
    for (idx r = 0; r < n; r++) lcp[r] = plcp[SA[r]];
    lcp[0] = 0;
    log("LCP array");
}

// ---------------------------------------------------------------------------
// tables
// ---------------------------------------------------------------------------
static const uint64_t EXACT = 1024;
static uint32_t bucket(uint64_t x) {  // exact up to 1024, then 1% geometric bins
    return x <= EXACT ? (uint32_t)x : (uint32_t)(EXACT + std::ceil(std::log((double)x / EXACT) / std::log(1.01)));
}
static uint64_t rep(uint32_t b) {  // representative value of a bin
    return b <= EXACT ? b : (uint64_t)std::llround(EXACT * std::pow(1.01, b - EXACT - 0.5));
}
static void add(std::vector<uint64_t>& v, size_t i, uint64_t x = 1) {
    if (v.size() <= i) v.resize(i + 1, 0);
    v[i] += x;
}
static void addf(std::vector<double>& v, size_t i, double x) {
    if (v.size() <= i) v.resize(i + 1, 0);
    v[i] += x;
}

static const uint32_t NB = 2600;  // > bucket(2^32)

// The P table of one radius r: cells (bin of c, bin of m); per cell the positions and, for
// e = 0..r, sum[e] = the sum over its positions of N_e, the number of copies (itself included)
// whose context agrees with its own on at least e bases on either side. N_0 = c, N_r = m, and
// sum[e] - sum[e + 1] counts the (position, copy) pairs whose contexts agree on exactly e bases.
struct Radius {
    int r;
    std::vector<int32_t> id;              // bin(c) * NB + bin(m) -> cell
    std::vector<uint32_t> cell;           // per k-mer position of G (c >= 2)
    std::vector<std::pair<uint32_t, uint32_t>> bins;
    std::vector<uint64_t> positions, sum; // per cell; sum: r + 1 per cell
};

struct Tables {
    std::vector<uint64_t> S, E;        // S by exact count, E by bin
    std::vector<double> F, C;          // families by exact count, copies by bin
    std::vector<Radius> P;
    std::vector<std::pair<int, std::vector<uint64_t>>> K;  // k-sweep: k' -> positions by exact count
    uint64_t positions = 0, repeated = 0;
    double distinct = 0, families = 0, copies = 0;
};

// Left extensions, from the suffix array of G reversed: the k-mer at G position i starts at
// g - i - k there, and the bases its copies share to its right are those they share to the left
// in G. Le[i] = UINT32_MAX where the k-mer occurs once (or is not a k-mer).
static std::vector<uint32_t> left_extensions(uint64_t g, int k, const std::vector<idx>& SA,
                                             const std::vector<uint32_t>& lcp) {
    std::vector<uint32_t> Le(g, UINT32_MAX);
    uint64_t a = 0;
    uint32_t rmin = UINT32_MAX;
    for (uint64_t r = 1; r <= g; r++) {
        if (r < g && lcp[r] >= (uint32_t)k) {
            rmin = std::min(rmin, lcp[r]);
            continue;
        }
        if (r - a >= 2)
            for (uint64_t q = a; q < r; q++) Le[g - (uint64_t)SA[q] - k] = rmin - k;
        a = r, rmin = UINT32_MAX;
    }
    return Le;
}

static void analyse(const std::vector<uint8_t>& T, int k, const std::vector<int>& radii, const std::vector<int>& ks,
                    const std::vector<idx>& SA, const std::vector<uint32_t>& lcp, std::vector<uint32_t>& Le, Tables& t) {
    const uint64_t g = T.size();
    std::vector<uint32_t> c(g), R(g, 0);  // per k-mer position of G
    uint64_t run = 0;
    for (uint64_t i = g; i-- > 0;) {  // valid k-mers: count 1 until an interval says otherwise
        run = T[i] ? run + 1 : 0;
        c[i] = run >= (uint64_t)k;
    }
    // k-mer intervals: maximal runs of suffixes with lcp >= k
    uint64_t a = 0;
    uint32_t rmin = UINT32_MAX;
    for (uint64_t r = 1; r <= g; r++) {
        if (r < g && lcp[r] >= (uint32_t)k) {
            rmin = std::min(rmin, lcp[r]);
            continue;
        }
        if (r - a >= 2) {
            uint32_t cnt = (uint32_t)std::min<uint64_t>(r - a, UINT32_MAX);
            for (uint64_t q = a; q < r; q++) c[SA[q]] = cnt, R[SA[q]] = rmin;
        }
        a = r, rmin = UINT32_MAX;
    }
    for (uint64_t i = 0; i < g; i++) {
        if (!c[i]) continue;
        t.positions++, add(t.S, c[i]), t.distinct += 1.0 / c[i];
        if (c[i] < 2) continue;
        t.repeated++;
        uint64_t extent = (uint64_t)Le[i] + R[i];
        double copy = 1.0 / (extent - k + 1);
        add(t.E, bucket(extent)), addf(t.C, bucket(extent), copy), addf(t.F, c[i], copy / c[i]);
        t.copies += copy, t.families += copy / c[i];
    }
    log("copy numbers, extents and copies");
    R = std::vector<uint32_t>(), Le = std::vector<uint32_t>();  // free
    // context counts: N_e(i) = occurrences of the (k + 2e)-mer centred on the k-mer at i, i.e. of
    // the suffix i - e cut at k + 2e bases (the size of its LCP interval); e from the largest
    // radius down to 0 (N_0 = c), so the cells of a radius are known before its smaller e
    int e_max = *std::max_element(radii.begin(), radii.end());
    for (int rad : radii) {
        Radius P;
        P.r = rad, P.id.assign((size_t)NB * NB, -1);
        t.P.push_back(std::move(P));
    }
    auto visit = [&](uint64_t i, int e, uint64_t N) {  // k-mer position i (c >= 2), its N_e
        for (Radius& P : t.P) {
            if (P.r < e) continue;
            if (P.r == e) {  // m = N_r: the position's cell
                size_t key = (size_t)bucket(c[i]) * NB + bucket(N);
                if (P.id[key] < 0) {
                    P.id[key] = (int32_t)P.bins.size();
                    P.bins.emplace_back(bucket(c[i]), bucket(N));
                    P.positions.push_back(0), P.sum.resize(P.sum.size() + P.r + 1, 0);
                }
                if (P.cell.empty()) P.cell.assign(g, 0);
                P.cell[i] = P.id[key], P.positions[P.id[key]]++;
            }
            P.sum[(size_t)P.cell[i] * (P.r + 1) + e] += N;
        }
    };
    for (int e = e_max; e >= 0; e--) {
        if (e == 0) {
            for (uint64_t i = 0; i < g; i++)
                if (c[i] >= 2) visit(i, 0, c[i]);
            break;
        }
        uint32_t ell = k + 2 * e;
        a = 0;
        for (uint64_t r = 1; r <= g; r++) {
            if (r < g && lcp[r] >= ell) continue;
            for (uint64_t q = a; q < r; q++) {
                uint64_t i = (uint64_t)SA[q] + e;
                if (i < g && c[i] >= 2) visit(i, e, r - a);
            }
            a = r;
        }
        for (uint64_t i = 0; i < std::min<uint64_t>(e, g); i++)  // no suffix starts e bases before them
            if (c[i] >= 2) visit(i, e, 1);
        if (std::find(radii.begin(), radii.end(), e) != radii.end()) log("context counts, radius " + std::to_string(e));
    }
    for (Radius& P : t.P) P.cell = std::vector<uint32_t>(), P.id = std::vector<int32_t>();
    // k-sweep: the spectrum of each k' from the intervals of suffixes sharing >= k' bases (all of them
    // k'-mers); the k'-mers counted once are the valid positions not in such an interval
    for (int kk : ks) {
        std::vector<uint64_t> S(2, 0);
        uint64_t valid = 0, in_rep = 0, len = 0;
        for (uint64_t i = 0; i <= g; i++) {  // N-free stretches of len bases hold len - k' + 1 k'-mers
            if (i < g && T[i]) {
                len++;
            } else {
                if (len >= (uint64_t)kk) valid += len - kk + 1;
                len = 0;
            }
        }
        uint64_t a = 0;
        for (uint64_t r = 1; r <= g; r++) {
            if (r < g && lcp[r] >= (uint32_t)kk) continue;
            if (r - a >= 2) add(S, r - a, r - a), in_rep += r - a;
            a = r;
        }
        S[1] = valid - in_rep;
        t.K.emplace_back(kk, std::move(S));
    }
    log("k-sweep spectra");
}

static void write(const std::string& out, const Genome& G, int k, const char* radii, const Tables& t) {
    FILE* f = fopen(out.c_str(), "w");
    fprintf(f, "#k\t%d\n#radii\t%s\n#records\t%ld\n#bases\t%llu\n", k, radii, G.records, (unsigned long long)G.acgt);
    fprintf(f, "#positions\t%llu\n#repeated\t%llu\n#distinct\t%.0f\n#c_max\t%zu\n#copies\t%.1f\n#families\t%.1f\n",
            (unsigned long long)t.positions, (unsigned long long)t.repeated, t.distinct, t.S.size() - 1, t.copies,
            t.families);
    fprintf(f, "#seconds\t%.1f\n#peak_gb\t%.2f\n", now() - t_start, peak_gb());
    for (size_t c = 1; c < t.S.size(); c++)
        if (t.S[c]) fprintf(f, "S\t%zu\t%llu\n", c, (unsigned long long)t.S[c]);
    for (auto& kv : t.K)
        for (size_t c = 1; c < kv.second.size(); c++)
            if (kv.second[c]) fprintf(f, "K\t%d\t%zu\t%llu\n", kv.first, c, (unsigned long long)kv.second[c]);
    for (const Radius& P : t.P) {
        std::vector<size_t> order(P.bins.size());
        for (size_t j = 0; j < order.size(); j++) order[j] = j;
        std::sort(order.begin(), order.end(), [&](size_t x, size_t y) { return P.bins[x] < P.bins[y]; });
        for (size_t j : order) {
            fprintf(f, "P\t%d\t%llu\t%llu\t%llu", P.r, (unsigned long long)rep(P.bins[j].first),
                    (unsigned long long)rep(P.bins[j].second), (unsigned long long)P.positions[j]);
            const uint64_t* sum = &P.sum[j * (P.r + 1)];
            for (int e = 0; e < P.r; e++) fprintf(f, "\t%llu", (unsigned long long)(sum[e] - sum[e + 1]));
            fprintf(f, "\n");
        }
    }
    for (size_t b = 0; b < t.E.size(); b++)
        if (t.E[b]) fprintf(f, "E\t%llu\t%llu\n", (unsigned long long)rep(b), (unsigned long long)t.E[b]);
    for (size_t c = 2; c < t.F.size(); c++)
        if (t.F[c] > 0) fprintf(f, "F\t%zu\t%.6g\n", c, t.F[c]);
    for (size_t b = 0; b < t.C.size(); b++)
        if (t.C[b] > 0) fprintf(f, "C\t%llu\t%.6g\n", (unsigned long long)rep(b), t.C[b]);
    fclose(f);
}

static std::vector<int> ints(const std::string& s) {  // "7,10" -> {7, 10}
    std::vector<int> v;
    std::stringstream ss(s);
    for (std::string x; std::getline(ss, x, ',');) v.push_back(std::stoi(x));
    return v;
}

int main(int argc, char** argv) {
    if (argc != 5) return fputs("usage: copies genome.fa.gz k radii ks\n", stderr), 1;
    t_start = now();
    std::string path = argv[1], out = path.substr(0, path.size() - strlen(".fa.gz")) + ".copies.tsv";
    int k = std::stoi(argv[2]);
    std::vector<int> radii = ints(argv[3]), ks = ints(argv[4]);

    Genome G = read_fasta(path.c_str());
    const std::vector<uint8_t>& T = G.text;
    log("read " + path + ": " + std::to_string(G.records) + " records, " + std::to_string(G.acgt) + " ACGT");

    std::vector<uint32_t> Le;
    {  // left extensions: suffix array of G reversed
        std::vector<uint8_t> Tr(T.rbegin(), T.rend());
        std::vector<idx> SA;
        std::vector<uint32_t> lcp;
        suffix_lcp(Tr, SA, lcp);
        Le = left_extensions(Tr.size(), k, SA, lcp);
        log("left extensions (reversed genome)");
    }
    Tables t;
    {
        std::vector<idx> SA;
        std::vector<uint32_t> lcp;
        suffix_lcp(T, SA, lcp);
        analyse(T, k, radii, ks, SA, lcp, Le, t);
    }
    write(out, G, k, argv[3], t);
    log("wrote " + out);
    return 0;
}
