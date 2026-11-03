// ue: uniqueness entropy of k-mer sketches.
//
// Reads FASTA files, passes every record to each sketch (sketches.cpp),
// computes the statistics and the Theorem 1 predictions (theory.cpp) and
// writes one JSON file. src/plot.py turns it into tables and figures.
//
// build: make -C src ue, or c++ -O3 -std=c++17 -ffp-contract=off -o src/ue src/main.cpp src/sketches.cpp src/theory.cpp -lz
//        (-ffp-contract=off: no fused multiply-add, so floating point matches numpy)
#include <zlib.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "sketches.h"
#include "theory.h"

static const char* USAGE = R"(usage: ue <command> [genome.fa ...] [--option value ...]

commands
  panel           every scheme at matched density: H, H_true, Delta, U_N, U_D,
                  E-hits, spurious mass, read anchoring, and Theorem 1 predictions
                  when <genome>.truth.tsv exists and --calibration is given
  closed-forms    Zipf closed forms vs the simulated repeat entropy (needs truth)
  estimate-alpha  alpha from a FracMinHash sketch (needs truth for the true alpha)
  real            real genomes: spectrum, panel, alpha from a sketch, k-sweep, and predictions
                  of H per scheme when <genome>.copies.tsv exists (written by src/copies)
  sketch          print the seeds of one scheme (--scheme) as TSV
  settings        print k, w, density, the mod-minimizer window and the context radii of
                  the schemes (the radii src/copies needs)

Each FASTA is one genome; its records are sketched separately and every
non-ACGT base splits them. Seeds are read on the forward strand only: a
k-mer and its reverse complement are different seeds. Settings come from
the command line, else from the genome's truth file (written by
src/simulate.py), else the defaults.

options [defaults]
  --out FILE         JSON output (stdout if omitted)
  --k INT            [21]
  --w INT            window; matched density 2/(w+1) [11]
  --density FLOAT    target density instead of 2/(w+1)
  --modmin-w INT     mod-minimizer window [searched for the density]
  --schemes A,B,...  [all,random-min,lex-min,weighted-min,mod-min,open-sync,
                     closed-sync,fracminhash]; also modmer, randstrobe, oracle
  --hash-seeds A,B   hash seeds [0,1,2; real: 0,1]
  --reads INT        error-free reads for anchoring statistics [0; real: 10000]
  --read-len INT     [1000; real: 150]
  --calibration FILE C and b of the context-dependent schemes, written by src/calibrate
                     (panel predictions need it)
  --no-predict       skip predictions
  --max-bases INT    read at most this many bases per FASTA
  --quiet            no progress messages
  estimate-alpha:  --sketch-density [0.05]
  real:            --ks [15,17,19,21,25,31] (--ks none to skip) --c-min [2] --sketch-density [0.05]
                   --no-tandem: predictions without the tandem term (periodic k-mers selected at rate d)
  sketch:          --scheme [random-min] --hash-seeds [0]
)";

static const char* PANEL = "all,random-min,lex-min,weighted-min,mod-min,open-sync,closed-sync,fracminhash";
static const double NaN = std::nan("");

// ---------------------------------------------------------------------------
// command line: --key value, --flag; lists are comma-separated
// ---------------------------------------------------------------------------
static std::map<std::string, std::string> opt;
static std::vector<std::string> files;
static const std::set<std::string> FLAGS = {"no-predict", "no-tandem", "quiet", "help"};
static const std::set<std::string> VALUES = {"out", "k", "w", "density", "modmin-w", "schemes", "hash-seeds",
                                             "reads", "read-len", "max-bases", "calibration",
                                             "sketch-density", "ks", "c-min", "scheme"};

static std::vector<std::string> split(const std::string& s, char sep = ',') {
    std::vector<std::string> out;
    std::stringstream ss(s);
    for (std::string x; std::getline(ss, x, sep);)
        if (!x.empty()) out.push_back(x);
    return out;
}

static void log(const char* fmt, ...) __attribute__((format(printf, 1, 2)));
static void log(const char* fmt, ...) {
    if (opt.count("quiet")) return;
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fputc('\n', stderr);
}

static double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// ---------------------------------------------------------------------------
// genomes: FASTA records plus the truth file of simulated genomes
// ---------------------------------------------------------------------------
struct Copy { long family, copies, length, start; };

struct Genome {
    std::string name;                     // file name without directory and extensions
    std::vector<std::string> names, seq;  // records
    std::vector<uint64_t> offset;         // start of each record in genome coordinates
    std::vector<std::pair<std::string, std::string>> truth;  // '#key value' lines of the truth file
    std::vector<Copy> copies;             // planted repeat copies (in the first record)
};

static std::string strip_ext(std::string p) {
    for (const char* e : {".gz", ".fa", ".fasta", ".fna"})
        if (p.size() > strlen(e) && p.compare(p.size() - strlen(e), strlen(e), e) == 0) p.resize(p.size() - strlen(e));
    return p;
}

static Genome load(const std::string& path) {
    Genome g;
    g.name = strip_ext(path.substr(path.find_last_of('/') + 1));
    gzFile f = gzopen(path.c_str(), "rb");  // plain or gzip
    if (!f) throw std::runtime_error("cannot open " + path);
    uint64_t max_bases = opt.count("max-bases") ? std::stoull(opt["max-bases"]) : 0, total = 0;
    static char buf[1 << 16];
    bool header = false, header_done = false, line_start = true;
    int n;
    while ((n = gzread(f, buf, sizeof buf)) > 0) {
        for (int i = 0; i < n; i++) {
            char c = buf[i];
            if (c == '\n') {
                header = header_done = false, line_start = true;
            } else if (line_start && c == '>') {
                header = true, line_start = false;
                g.names.emplace_back(), g.seq.emplace_back();
            } else if (header) {
                if ((c == ' ' || c == '\t') && !g.names.back().empty()) header_done = true;  // name: first word
                if (!header_done && c > ' ') g.names.back() += c;
            } else if (!g.seq.empty() && c > ' ' && (!max_bases || total < max_bases)) {
                g.seq.back() += c, total++, line_start = false;
            } else {
                line_start = false;
            }
        }
    }
    gzclose(f);
    while (!g.seq.empty() && g.seq.back().empty()) g.seq.pop_back(), g.names.pop_back();
    if (g.seq.empty()) throw std::runtime_error("no sequence in " + path);
    g.offset.push_back(0);
    for (auto& s : g.seq) g.offset.push_back(g.offset.back() + s.size());

    std::ifstream t(strip_ext(path) + ".truth.tsv");
    for (std::string line; std::getline(t, line);) {
        if (line.empty() || line.rfind("family", 0) == 0) continue;
        if (line[0] == '#') {
            size_t tab = line.find('\t');
            g.truth.emplace_back(line.substr(1, tab - 1), line.substr(tab + 1));
        } else {
            Copy c;
            std::istringstream(line) >> c.family >> c.copies >> c.length >> c.start;
            g.copies.push_back(c);
        }
    }
    log("[%s] %zu record(s), %llu bases%s", g.name.c_str(), g.seq.size(), (unsigned long long)g.offset.back(),
        g.truth.empty() ? "" : " + truth");
    return g;
}

static bool has_truth(const Genome& g) { return !g.truth.empty(); }

// a setting: the command line, else the truth file, else the default
static std::string setting(const Genome* g, const std::string& key, const std::string& def) {
    if (opt.count(key)) return opt[key];
    std::string tkey = key;
    std::replace(tkey.begin(), tkey.end(), '-', '_');
    if (g)
        for (auto& kv : g->truth)
            if (kv.first == tkey) return kv.second;
    return def;
}

// ---------------------------------------------------------------------------
// JSON output: a row is a list of (key, JSON value)
// ---------------------------------------------------------------------------
typedef std::vector<std::pair<std::string, std::string>> Row;

static std::string num(double x) {
    if (std::isnan(x)) return "NaN";  // read back by Python's json module
    if (std::isinf(x)) return x > 0 ? "Infinity" : "-Infinity";
    char b[32];
    for (int p = 15; p <= 17; p++) {  // shortest that round-trips
        snprintf(b, sizeof b, "%.*g", p, x);
        if (strtod(b, nullptr) == x) break;
    }
    return b;
}

static std::string quote(const std::string& s) {
    std::string o = "\"";
    for (char c : s) o += (c == '"' || c == '\\') ? std::string("\\") + c : std::string(1, c);
    return o + "\"";
}

static std::string value(const std::string& s) {  // truth values: numbers stay numbers
    if (s == "None") return "null";
    char* end;
    strtod(s.c_str(), &end);
    return (!s.empty() && *end == 0) ? s : quote(s);
}

static void put(Row& r, const std::string& key, const std::string& json) {
    for (auto& kv : r)
        if (kv.first == key) { kv.second = json; return; }
    r.emplace_back(key, json);
}
static void put(Row& r, const std::string& key, double x) { put(r, key, num(x)); }
static void put_str(Row& r, const std::string& key, const std::string& s) { put(r, key, quote(s)); }

static std::string get_str(const Row& r, const std::string& key) {
    for (auto& kv : r)
        if (kv.first == key) return kv.second.substr(1, kv.second.size() - 2);
    return "";
}

static double get(const Row& r, const std::string& key) {
    for (auto& kv : r)
        if (kv.first == key) return strtod(kv.second.c_str(), nullptr);
    return NaN;
}

static std::string array(const std::vector<double>& v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); i++) s += (i ? ", " : "") + num(v[i]);
    return s + "]";
}

// {"command": ..., "<table>": [rows], ..., <extra>}
static void write_json(const std::string& command, const std::vector<std::pair<std::string, std::vector<Row>>>& tables,
                       const std::string& extra = "") {
    std::string out = opt.count("out") ? opt["out"] : "";
    if (!out.empty() && std::filesystem::path(out).has_parent_path())
        std::filesystem::create_directories(std::filesystem::path(out).parent_path());
    FILE* f = out.empty() ? stdout : fopen((out + ".tmp").c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + out);
    fprintf(f, "{\"command\": %s", quote(command).c_str());
    for (auto& t : tables) {
        fprintf(f, ",\n\"%s\": [", t.first.c_str());
        for (size_t i = 0; i < t.second.size(); i++) {
            fprintf(f, "%s\n {", i ? "," : "");
            for (size_t j = 0; j < t.second[i].size(); j++)
                fprintf(f, "%s%s: %s", j ? ", " : "", quote(t.second[i][j].first).c_str(), t.second[i][j].second.c_str());
            fprintf(f, "}");
        }
        fprintf(f, "\n]");
    }
    fprintf(f, "%s}\n", extra.c_str());
    if (f != stdout) {
        fclose(f);
        rename((out + ".tmp").c_str(), out.c_str());
        log("wrote %s", out.c_str());
    }
}

// ---------------------------------------------------------------------------
// the panel: every scheme at matched density d. Settings, select_seeds(), seed_space(), radius()
// and the other scheme parameters are in sketches.cpp
// ---------------------------------------------------------------------------
static Settings settings(const Genome* g) {
    Settings s;
    s.k = std::stoi(setting(g, "k", "21"));
    s.w = std::stoi(setting(g, "w", "11"));
    s.d = std::stod(setting(g, "density", num(2.0 / (s.w + 1))));
    // the truth file's window was chosen for its k and w: ignore it when they are overridden
    bool overridden = opt.count("k") || opt.count("w") || opt.count("density");
    std::string mw = opt.count("modmin-w") ? opt["modmin-w"] : overridden ? "" : setting(g, "modmin-w", "");
    s.modmin_w = mw.empty() ? modmin_window(s.k, s.d) : std::stoi(mw);
    s.genome = nullptr;
    s.c_star = 0;
    s.tie = 0;
    return s;
}

// counts of every seed of a seed space over all records
static Counts count_space(const Genome& g, const std::string& key, const Settings& s, uint64_t hs) {
    std::vector<uint64_t> h;
    for (auto& rec : g.seq)
        for (const Seed& x : seed_space(key, rec, s, hs)) h.push_back(x.hash);
    return make_counts(std::move(h));
}

// positions per multiplicity: hist[c] = positions whose seed occurs c times
static std::vector<double> histogram(const Counts& counts) {
    const std::vector<uint64_t>& sorted = counts.sorted;
    std::vector<double> hist;
    for (size_t i = 0; i < sorted.size();) {
        size_t j = i;
        while (j < sorted.size() && sorted[j] == sorted[i]) j++;
        if (hist.size() <= j - i) hist.resize(j - i + 1, 0);
        hist[j - i] += j - i;
        i = j;
    }
    return hist;
}

// ---------------------------------------------------------------------------
// statistics (ue.metrics.measure); ms/mg: sketch and genome multiplicity of
// every selected seed; rep: seed lies in a planted repeat copy (optional)
// ---------------------------------------------------------------------------
static Row measure(const std::vector<uint32_t>& ms, const std::vector<uint32_t>& mg, const std::vector<uint8_t>* rep,
                   double n_valid, int k_eff, int r) {
    Row out;
    double N = ms.size();
    if (N == 0) {
        put(out, "N", 0), put(out, "d", 0);
        return out;
    }
    std::vector<double> hs, hg;
    double sing = 0, spur = 0;
    for (size_t i = 0; i < ms.size(); i++) {
        if (hs.size() <= ms[i]) hs.resize(ms[i] + 1, 0);
        if (hg.size() <= mg[i]) hg.resize(mg[i] + 1, 0);
        hs[ms[i]]++, hg[mg[i]]++;
        if (ms[i] == 1) sing++, spur += mg[i] > 1;
    }
    double H = 0, H_true = 0, distinct = 0, eh = 0, eht = 0;
    for (size_t c = 1; c < hs.size(); c++) H += hs[c] * std::log2(c), distinct += hs[c] / c, eh += hs[c] * c;
    for (size_t c = 1; c < hg.size(); c++) H_true += hg[c] * std::log2(c), eht += hg[c] * c;
    H /= N, H_true /= N;
    double d = N / n_valid;
    put(out, "N", N), put(out, "d", d), put(out, "k_eff", k_eff), put(out, "r", r);
    put(out, "H", H);                          // uniqueness entropy
    put(out, "H_true", H_true);                // with genome multiplicities
    put(out, "U_N", hs.size() > 1 ? hs[1] / N : 0);         // singleton fraction (positions)
    put(out, "U_D", hs.size() > 1 ? hs[1] / distinct : 0);  // singleton fraction (distinct seeds)
    put(out, "Ehits", eh / N), put(out, "Ehits_true", eht / N);
    put(out, "H_min", std::log2(hs.size() - 1.0));          // log2 of the largest multiplicity
    put(out, "spurious_mass", spur / N);
    put(out, "Delta", H_true - H);
    put(out, "singleton_precision", sing ? 1 - spur / sing : NaN);
    if (rep) {  // background and repeat parts
        double n[2] = {0, 0}, h[2] = {0, 0}, ht[2] = {0, 0};
        for (size_t i = 0; i < ms.size(); i++) {
            int t = (*rep)[i];
            n[t]++, h[t] += std::log2(ms[i]), ht[t] += std::log2(mg[i]);
        }
        const char* tag[2] = {"bg", "rep"};
        for (int t = 0; t < 2; t++) {
            put(out, std::string("frac_") + tag[t], n[t] / N);
            put(out, std::string("H_") + tag[t], n[t] ? h[t] / n[t] : NaN);
            put(out, std::string("H_true_") + tag[t], n[t] ? ht[t] / n[t] : NaN);
        }
    }
    return out;
}

// error-free reads sampled uniformly within records (ue.experiment.read_anchor_stats);
// pos: genome coordinates of the selected seeds, ascending
static void read_stats(Row& row, const Genome& g, const std::vector<uint64_t>& pos, const std::vector<uint32_t>& ms,
                       const std::vector<uint32_t>& mg, int k_eff, int span, int n_reads, int read_len, uint64_t seed) {
    long L = std::max(1, read_len - k_eff + 1);  // seed starts inside a read
    std::vector<uint64_t> cum = {0};
    for (auto& s : g.seq) cum.push_back(cum.back() + std::max(0L, (long)s.size() - span + 1 - L));
    if (cum.back() == 0) return;
    std::mt19937_64 rng(seed);
    double anchors = 0, hits = 0, tu_reads = 0, fu_reads = 0, none = 0;
    for (int i = 0; i < n_reads; i++) {
        uint64_t x = rng() % cum.back();
        size_t r = std::upper_bound(cum.begin(), cum.end(), x) - cum.begin() - 1;
        uint64_t a = g.offset[r] + (x - cum[r]);
        long tu = 0, fu = 0;
        for (size_t j = std::lower_bound(pos.begin(), pos.end(), a) - pos.begin(); j < pos.size() && pos[j] < a + L; j++) {
            anchors++, hits += ms[j];
            if (ms[j] == 1) (mg[j] == 1 ? tu : fu)++;
        }
        if (tu) tu_reads++;
        else if (fu) fu_reads++;
        else none++;
    }
    put(row, "unique_anchor_rate", tu_reads / n_reads);
    put(row, "false_unique_rate", fu_reads / n_reads);
    put(row, "no_unique_rate", none / n_reads);
    put(row, "anchors_per_read", anchors / n_reads);
    put(row, "hits_per_read", hits / n_reads);
    put(row, "hits_per_anchor", hits / std::max(anchors, 1.0));
}

// P(selected | planted copy number c); rates[1] is the background
static std::vector<double> selection_function(const Genome& g, const std::vector<uint64_t>& pos,
                                              const std::vector<uint32_t>& planted, int k, double& overall) {
    std::vector<uint8_t> sel(g.offset.back(), 0);
    for (uint64_t p : pos) sel[p] = 1;
    std::vector<double> tot, hit;
    double all = 0, all_hit = 0;
    for (size_t r = 0; r < g.seq.size(); r++) {
        long run = 0;  // valid bases ending here
        for (size_t i = 0; i < g.seq[r].size(); i++) {
            run = strchr("ACGTacgt", g.seq[r][i]) && g.seq[r][i] ? run + 1 : 0;
            if (run < k) continue;
            size_t p = g.offset[r] + i + 1 - k;
            uint32_t c = p < planted.size() ? planted[p] : 0;
            if (tot.size() <= c) tot.resize(c + 1, 0), hit.resize(c + 1, 0);
            tot[c]++, hit[c] += sel[p], all++, all_hit += sel[p];
        }
    }
    overall = all_hit / all;
    std::vector<double> rate(std::max<size_t>(tot.size(), 2), NaN);
    rate[1] = tot.size() && tot[0] ? hit[0] / tot[0] : 1;
    for (size_t c = 2; c < tot.size(); c++)
        if (tot[c]) rate[c] = hit[c] / tot[c];
    return rate;
}

// ---------------------------------------------------------------------------
// predictions (ue.experiment.predictions)
// ---------------------------------------------------------------------------
// C, b, b by copy length and the agreement of chance matches of the context-dependent schemes, read
// from the file src/calibrate writes (--calibration): one row per scheme, k, w and theta
struct Calib {
    std::string scheme;
    int k, w, r;
    double theta, C, b, coll;
    std::vector<double> L, b_L;  // b of copies of length L; copies longer than the last L have b
};
static std::vector<Calib> calibration;

static void load_calibration(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::string line;
    std::getline(f, line);
    std::vector<std::string> head = split(line, '\t');
    auto field = [&head](const std::vector<std::string>& row, const std::string& name) {
        size_t i = std::find(head.begin(), head.end(), name) - head.begin();
        if (i >= row.size()) throw std::runtime_error("calibration file: no column " + name);
        return row[i];
    };
    while (std::getline(f, line)) {
        std::vector<std::string> row = split(line, '\t');
        if (row.empty()) continue;
        Calib c;
        c.scheme = field(row, "scheme");
        c.k = std::stoi(field(row, "k")), c.w = std::stoi(field(row, "w")), c.r = std::stoi(field(row, "r"));
        c.theta = std::stod(field(row, "theta")), c.C = std::stod(field(row, "C")), c.b = std::stod(field(row, "b"));
        c.coll = std::stod(field(row, "collision"));
        for (auto& x : split(field(row, "b_by_length"))) {  // L:b,L:b,...
            size_t colon = x.find(':');
            c.L.push_back(std::stod(x.substr(0, colon))), c.b_L.push_back(std::stod(x.substr(colon + 1)));
        }
        calibration.push_back(c);
    }
    log("read %zu calibration rows from %s", calibration.size(), path.c_str());
}

static const Calib& calibration_for(const std::string& key, const Settings& s, double theta) {
    for (const Calib& c : calibration)
        if (c.scheme == key && c.k == s.k && c.w == s.w && c.r == radius(key, s) && std::fabs(c.theta - theta) < 1e-9)
            return c;
    char msg[200];
    snprintf(msg, sizeof msg, "no calibration of %s at k=%d w=%d theta=%g: run src/calibrate with them and pass --calibration",
             key.c_str(), s.k, s.w, theta);
    throw std::runtime_error(msg);
}

static double interp(const std::vector<double>& xs, const std::vector<double>& ys, double x) {  // np.interp
    if (x <= xs.front()) return ys.front();
    if (x >= xs.back()) return ys.back();
    size_t j = std::upper_bound(xs.begin(), xs.end(), x) - xs.begin();
    return ys[j - 1] + (x - xs[j - 1]) / (xs[j] - xs[j - 1]) * (ys[j] - ys[j - 1]);
}

static void predictions(std::vector<Row>& rows, const Genome& g, const Settings& s,
                        const std::map<std::string, std::pair<std::vector<double>, double>>& d_of_c) {
    double theta = std::stod(setting(&g, "theta", "0"));
    std::map<long, std::pair<long, long>> fam_map;  // family -> (copies, length)
    for (auto& c : g.copies) fam_map[c.family] = {c.copies, c.length};
    Families fam;
    for (auto& f : fam_map) fam.push_back(f.second);
    double n = g.seq[0].size(), m = n - s.k + 1;
    for (Row& row : rows) {
        std::string key = get_str(row, "scheme");
        double d = get(row, "d");
        put(row, "H_pred", NaN), put(row, "H_pred_crude", NaN);
        if (key == "randstrobe" || key == "oracle") continue;
        Scheme sc;
        sc.k = s.k;
        if (radius(key, s) > 0) {  // context-dependent: C and b from the calibration file
            const Calib& c = calibration_for(key, s, theta);
            sc.r = c.r, sc.C = c.C, sc.b = c.b, sc.coll_agree = c.coll;
            sc.b_of_L = [c](double L) { return L <= c.L.back() ? interp(c.L, c.b_L, L) : c.b; };
            if (d_of_c.count(key)) {
                auto f = d_of_c.at(key);
                sc.d_of_c = [f](long cc) { return cc < (long)f.first.size() && !std::isnan(f.first[cc]) ? f.first[cc] : f.second; };
            }
            put(row, "C_cal", c.C), put(row, "b_cal", c.b);
            if (key == "random-min") {  // earlier form of Eq. (master) (no seed survival, no shared substitutions), crude C and b; not in the paper
                Scheme crude;
                crude.k = s.k, crude.r = sc.r, crude.C = crude_C(theta, s.w, d), crude.b = crude_b(d);
                put(row, "H_pred_crude", predict(crude, fam, m, n, theta, false, false, false).H);
            }
        }
        Prediction P = predict(sc, fam, m, n, theta);
        put(row, "H_pred", P.H), put(row, "H_genome_pred", P.H_genome), put(row, "Delta_ctx_pred", P.Delta_ctx);
        put(row, "Delta_cons_pred", P.Delta_cons), put(row, "Delta_order_pred", P.Delta_order);
    }
}

// ---------------------------------------------------------------------------
// predictions for a real genome: its repeat structure comes from src/copies (suffix array and
// LCP array), read from <genome>.copies.tsv; see the top of copies.cpp for the tables
// ---------------------------------------------------------------------------
struct CopyTables {
    std::string path;
    std::vector<std::pair<std::string, std::string>> header;    // '#key value' lines
    std::map<long, double> S, E, F, C;                          // spectrum, extents, families by c, copies by length
    struct Cell { double c, m, positions; std::vector<double> s; };  // s[e]: copies agreeing on exactly e bases
    std::map<int, std::vector<Cell>> P;                         // radius -> cells
    std::map<int, std::map<long, double>> K;                    // k-sweep: k -> (count -> positions)
    std::string get(const std::string& key) const {
        for (auto& kv : header)
            if (kv.first == key) return kv.second;
        return "";
    }
};

static bool load_copies(const std::string& fasta, CopyTables& t) {
    t.path = strip_ext(fasta) + ".copies.tsv";
    std::ifstream in(t.path);
    if (!in) return false;
    for (std::string line; std::getline(in, line);) {
        if (line.empty()) continue;
        if (line[0] == '#') {
            size_t tab = line.find('\t');
            t.header.emplace_back(line.substr(1, tab - 1), tab == std::string::npos ? "" : line.substr(tab + 1));
            continue;
        }
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "P") {
            int r;
            CopyTables::Cell cell;
            ss >> r >> cell.c >> cell.m >> cell.positions;
            for (double x; ss >> x;) cell.s.push_back(x);
            t.P[r].push_back(cell);
        } else if (tag == "K") {
            int kk;
            long c;
            double x;
            ss >> kk >> c >> x;
            t.K[kk][c] += x;
        } else if (tag == "S" || tag == "E" || tag == "F" || tag == "C") {
            long a;
            double x;
            ss >> a >> x;
            (tag == "S" ? t.S : tag == "E" ? t.E : tag == "F" ? t.F : t.C)[a] += x;
        }
    }
    return true;
}

static uint32_t bucket(uint64_t x) {  // the bins of copies.cpp: exact up to 1024, then 1% geometric
    return x <= 1024 ? (uint32_t)x : (uint32_t)(1024 + std::ceil(std::log(x / 1024.0) / std::log(1.01)));
}

static double weighted_median(const std::map<long, double>& h) {
    double total = 0, cum = 0;
    for (auto& x : h) total += x.second;
    for (auto& x : h)
        if ((cum += x.second) >= total / 2) return x.first;
    return NaN;
}

// ---------------------------------------------------------------------------
// tandem term: k-mers of minimal period p <= k/2 (inside tandem runs) are selected by leftmost-argmin
// schemes at the run rate 1/p, not at the density (Lemma: ties in tandem repeats). The genome enters
// through where its periodic k-mers sit (period, count c, distance to both ends of their run), each
// scheme through its run rates (tandem_rates), and Theorem 1 through d_S(c).
// ---------------------------------------------------------------------------
struct TandemTable {
    int dmax = 0, p_max = 0;
    std::map<uint64_t, double> T;         // key(bucket(c), p, dl, dr) -> positions
    std::vector<double> by_p, logc_by_p;  // periodic positions and their sum of log2 c, per period
    std::map<long, double> by_c;          // periodic positions per count c (exact)
    double positions = 0;                 // all k-mer positions
    static uint64_t key(uint32_t b, int p, int dl, int dr) {
        return ((uint64_t)b << 24) | ((uint64_t)p << 16) | ((uint64_t)dl << 8) | (uint64_t)dr;
    }
};

// dmax: the largest context radius of the panel (distances to the run ends are capped there)
static TandemTable tandem_table(const Genome& g, int k, const Counts& kmers, int dmax) {
    TandemTable t;
    t.dmax = dmax, t.p_max = k / 2;
    t.by_p.assign(t.p_max + 1, 0), t.logc_by_p.assign(t.p_max + 1, 0);
    t.positions = kmers.sorted.size();
    for (auto& rec : g.seq) {
        std::vector<PeriodicKmer> pk = periodic_kmers(rec, k, t.p_max);
        for (size_t lo = 0; lo < pk.size();) {  // a run: consecutive periodic k-mers of one period
            size_t hi = lo;
            while (hi + 1 < pk.size() && pk[hi + 1].start == pk[hi].start + 1 && pk[hi + 1].p == pk[lo].p) hi++;
            for (size_t j = lo; j <= hi; j++) {
                uint32_t c = count(kmers, pk[j].hash);
                int dl = std::min<int>(j - lo, dmax), dr = std::min<int>(hi - j, dmax);
                t.T[TandemTable::key(bucket(c), pk[j].p, dl, dr)]++;
                t.by_p[pk[j].p]++, t.logc_by_p[pk[j].p] += std::log2(c), t.by_c[c]++;
            }
            lo = hi + 1;
        }
    }
    return t;
}

// schemes whose selection in a run follows from their hash order. All k-mers, FracMinHash and modmers
// select every k-mer at rate d whatever its period; the lex and weighted minimizers use their measured
// selection function; randstrobes and the oracle are not predicted.
static bool tandem_scheme(const std::string& key) {
    return key == "random-min" || key == "mod-min" || key == "open-sync" || key == "closed-sync";
}

static std::map<std::string, TandemRates> tandem_cache;
static const TandemRates& tandem_calibration(const std::string& key, const Settings& s) {
    char id[160];
    snprintf(id, sizeof id, "%s %d %d %d %.17g", key.c_str(), s.k, s.w, s.modmin_w, s.d);
    auto it = tandem_cache.find(id);
    if (it != tandem_cache.end()) return it->second;
    double t0 = now();
    SketchHS sel = [key, s](const std::string& seq, uint64_t hs) { return select_seeds(key, seq, s, hs); };
    TandemRates t = tandem_rates(sel, radius(key, s), s.k, s.k / 2);
    log("    tandem rates %s in %.1fs: base=%.3f, interior rate / base: p=1 %.2f, p=2 %.2f, p=3 %.2f", key.c_str(),
        now() - t0, t.base, t.rel(1, t.r, t.r), t.rel(2, t.r, t.r), t.rel(3, t.r, t.r));
    return tandem_cache[id] = t;
}

// d_S(c) / d_S: aperiodic positions count 1, periodic ones their run rate relative to the density
static std::function<double(long)> tandem_tilt(const TandemTable& tt, const TandemRates& tr, const std::map<long, double>& S) {
    std::map<uint32_t, double> all, per, sel;  // bucket of c -> positions, periodic positions, their weight
    for (auto& x : S) all[bucket(x.first)] += x.second;
    for (auto& x : tt.T) {
        uint32_t b = (uint32_t)(x.first >> 24);
        int p = (x.first >> 16) & 0xff, dl = (x.first >> 8) & 0xff, dr = x.first & 0xff;
        per[b] += x.second, sel[b] += x.second * tr.rel(p, dl, dr);
    }
    std::map<uint32_t, double> tilt;
    for (auto& a : all)
        if (a.second > 0) tilt[a.first] = (a.second - per[a.first] + sel[a.first]) / a.second;
    return [tilt](long c) {
        auto it = tilt.find(bucket(c));
        return it != tilt.end() ? it->second : 1.0;
    };
}

// One row per genome: what copies found, and whether its spectrum is ue's.
static Row copies_summary(const std::string& name, const CopyTables& cp, const std::vector<double>& khist) {
    Row r;
    put_str(r, "genome", name);
    for (const char* key : {"k", "records", "bases", "positions", "repeated", "distinct", "c_max", "copies", "families",
                            "seconds", "peak_gb"})
        put(r, key, std::stod(cp.get(key).empty() ? "nan" : cp.get(key)));
    put_str(r, "radii", cp.get("radii"));
    put(r, "rho_repeated", get(r, "repeated") / get(r, "positions"));
    double L = 0, n = 0, Lc = 0, nc = 0;
    for (auto& e : cp.E) L += e.first * e.second, n += e.second;
    for (auto& e : cp.C) Lc += e.first * e.second, nc += e.second;
    put(r, "extent_mean", L / n), put(r, "extent_median", weighted_median(cp.E));       // position-weighted
    put(r, "copy_length_mean", Lc / nc), put(r, "copy_length_median", weighted_median(cp.C));  // per copy
    for (auto& pr : cp.P) {  // boundary fraction: repeated positions whose copies do not all share the context
        double rep = 0, bd = 0;
        for (auto& x : pr.second) rep += x.positions, bd += x.c != x.m ? x.positions : 0;
        put(r, "boundary_fraction_r" + std::to_string(pr.first), bd / rep);
    }
    double diff = 0;  // the spectrum of copies (suffix array) against ue's (hashes)
    for (size_t c = 1; c < std::max(khist.size(), (size_t)cp.S.rbegin()->first + 1); c++) {
        double a = c < khist.size() ? khist[c] : 0, b = cp.S.count(c) ? cp.S.at(c) : 0;
        diff += std::fabs(a - b);
    }
    put(r, "spectrum_mismatch", diff);  // positions counted differently; 0 when both agree
    return r;
}

// Theorem 1 under the planted-Zipf model with tandem runs: the aperiodic k-mers follow a Zipf law fitted
// to their positions (alpha, rho, the extents as copy lengths), the periodic ones are taken as measured:
// share mu_p and log2 counts, selected at the run rates of tr (null: rate d). Periodic and aperiodic
// k-mers are different seeds, so H is the mixture of the two parts, weighted by the selected shares.
static double zipf_with_tandem(const Scheme& sc, const TandemTable& tt, const TandemRates* tr, const CopyTables& cp,
                               long c_min, double N) {
    std::vector<double> hist;  // aperiodic positions per count
    for (auto& x : cp.S) {
        if (hist.size() <= (size_t)x.first) hist.resize(x.first + 1, 0);
        hist[x.first] += x.second;
    }
    for (auto& x : tt.by_c)
        if ((size_t)x.first < hist.size()) hist[x.first] = std::max(0.0, hist[x.first] - x.second);
    double n_ap = 0, rep = 0;
    long cmax = 0;
    for (size_t c = 1; c < hist.size(); c++) {
        n_ap += hist[c], rep += c >= 2 ? hist[c] : 0;
        if (hist[c] > 0) cmax = c;
    }
    double a = fit_zipf_mle(hist, std::max(2L, c_min), true, 0);
    double H_ap = std::isfinite(a) && cmax >= 2 ? predict_zipf(sc, a, cmax, cp.E, rep / n_ap, N) : 0;
    double w = 0, h = 0;  // periodic part: selection weight relative to aperiodic positions, and log2 counts
    for (auto& x : tt.T) {
        uint32_t b = (uint32_t)(x.first >> 24);
        int p = (x.first >> 16) & 0xff, dl = (x.first >> 8) & 0xff, dr = x.first & 0xff;
        double rel = tr ? tr->rel(p, dl, dr) : 1;
        // representative count of the bin (exact up to 1024)
        double c = b <= 1024 ? b : 1024 * std::pow(1.01, b - 1024);
        w += x.second * rel, h += x.second * rel * std::log2(c);
    }
    double tau = w / (n_ap + w);  // share of selected positions in tandem runs
    return (1 - tau) * H_ap + (w > 0 ? tau * h / w : 0);
}

// Predicted H of each scheme from the genome's repeat structure, theta = 0:
//   H_pred       per position: of the c copies of its k-mer, the m whose whole context of radius r
//                is identical are selected with it; a copy whose context agrees on exactly e < r
//                bases on either side is selected with probability q(e) (context_kernel); with q
//                the mean of q(e) over these copies, E[log2(m + Bin(c - m, q))]
//   H_pred_thm   Theorem 1 with the genome's own spectrum: a position whose copies all share the
//                context (m = c) counts log2 c, the others (the boundary fraction) Phi(c - 1, b)
//   H_pred_zipf  Theorem 1 under the Zipf law fitted to the k-mer spectrum (alpha, c_max), with the
//                repeated fraction rho and the extents as copy lengths, beta(L) = 2r / (L - k + 1)
// theta = 0, so C = 1; b, q(e) and the collision agreement are calibrated on random sequence pairs
// (real_calibration). Positions are weighted by the selection rate d(c) of Table 1: measured for the lex and
// weighted minimizers, constant for the others. H_pred_measured_d is H_pred with every scheme's
// measured d(c), which separates false uniqueness (Delta) from selection that is not hash-random.
// Randstrobes and the oracle are not predicted.
// What a real-genome prediction needs from a scheme of radius r > 0, measured once per scheme on
// random sequences. C is not among them: at theta = 0 copies with an identical context are always
// selected together (C = 1).
struct RealCalib {
    double b, coll;        // edge agreement b (Theorem 1), agreement of chance collisions
    std::vector<double> q;  // context kernel q(e), e = 0..r-1
};
static std::map<std::string, RealCalib> real_calib_cache;

static const RealCalib& real_calibration(const std::string& key, const Settings& s) {
    char id[160];
    snprintf(id, sizeof id, "%s %d %d %d", key.c_str(), s.k, s.w, s.modmin_w);
    auto it = real_calib_cache.find(id);
    if (it != real_calib_cache.end()) return it->second;
    double t0 = now();
    Sketch sel = calibration_sketch(key, s), sp = [key, s](const std::string& seq) { return seed_space(key, seq, s, 0); };
    int r = radius(key, s);
    RealCalib c;
    c.b = calibrate(sel, sp, r, span(key, s), s.k, 0, std::max(1500, 4 * s.w), 12).b;
    c.coll = collision_agreement(sel, r, s.k, 3000);
    c.q = context_kernel(sel, r, s.k);
    log("    calibrated %s in %.1fs: b=%.3f q(0)=%.3f q(%d)=%.3f", key.c_str(), now() - t0, c.b, c.q[0], r - 1, c.q.back());
    return real_calib_cache[id] = c;
}

static void real_predictions(std::vector<Row>& rows, const Settings& s, const CopyTables& cp, const Row& spec,
                             const std::map<std::string, std::pair<std::vector<double>, double>>& d_of_c,
                             std::vector<Row>& kernels, const TandemTable* tt) {
    // Zipf fit: the position-weighted alpha is the one H depends on; the distinct-seed fit is kept for comparison
    double alpha = get(spec, "alpha_pos_trunc"), alpha_d = get(spec, "alpha_mle_trunc");
    double rho = get(spec, "rho_repeated"), N = get(spec, "positions");
    long c_max = (long)get(spec, "c_max"), c_min = std::stol(setting(nullptr, "c-min", "2"));
    auto rate_of = [&](const std::string& key) -> std::function<double(long)> {  // measured, pooled over bins
        if (!d_of_c.count(key)) return [](long) { return 1.0; };
        const std::vector<double>& rate = d_of_c.at(key).first;
        std::map<uint32_t, std::pair<double, double>> pool;  // bin -> (selected, positions)
        for (auto& x : cp.S)
            if (x.first < (long)rate.size() && !std::isnan(rate[x.first]))
                pool[bucket(x.first)].first += rate[x.first] * x.second, pool[bucket(x.first)].second += x.second;
        double overall = d_of_c.at(key).second;
        return [pool, overall](long c) {
            auto it = pool.find(bucket(c));
            return it != pool.end() && it->second.second > 0 ? it->second.first / it->second.second : overall;
        };
    };
    for (Row& row : rows) {
        std::string key = get_str(row, "scheme");
        for (const char* col : {"H_pred", "H_pred_thm", "H_pred_zipf", "H_pred_zipf_distinct", "H_pred_measured_d",
                                "H_pred_no_tandem", "H_true_pred", "H_true_pred_no_tandem", "Delta_pred", "beta_emp",
                                "beta_zipf", "agree_pred"})
            put(row, col, NaN);
        put(row, "alpha_fit", alpha), put(row, "alpha_fit_distinct", alpha_d), put(row, "rho_fit", rho), put(row, "c_max_fit", c_max);
        if (key == "randstrobe" || key == "oracle") continue;
        Scheme sc;
        sc.k = s.k, sc.r = radius(key, s);
        if (sc.r > 0 && !cp.P.count(sc.r)) {
            log("    %s: no context counts for radius %d in %s (copies radii)", key.c_str(), sc.r, cp.path.c_str());
            continue;
        }
        std::vector<double> q;  // context kernel q(e), e < r
        if (sc.r > 0) {
            const RealCalib& c = real_calibration(key, s);
            sc.C = 1, sc.b = c.b, sc.coll_agree = c.coll, q = c.q;
            put(row, "b_cal", c.b);
            bool listed = false;
            for (auto& kr : kernels) listed |= get_str(kr, "scheme") == key;
            if (!listed) {
                Row kr;
                put_str(kr, "scheme", key), put(kr, "r", sc.r), put(kr, "q", array(q)), put(kr, "b_cal", c.b);
                kernels.push_back(kr);
            }
        }
        // selection function d_S(c): measured for the lex and weighted minimizers (Table 1); from the tandem
        // term for schemes that break ties in runs; constant otherwise
        bool tilt = key == "lex-min" || key == "weighted-min", tandem = !tilt && tt && tandem_scheme(key);
        std::function<double(long)> one = [](long) { return 1.0; }, d_meas = rate_of(key);
        std::function<double(long)> d = tilt ? d_meas : tandem ? tandem_tilt(*tt, tandem_calibration(key, s), cp.S) : one;
        if (tilt || tandem) sc.d_of_c = d;
        // sums over positions weighted by d_S(c): [0] the model, [1] measured, [2] constant (no tandem term)
        std::function<double(long)> dj[3] = {d, d_meas, tilt ? d : one};
        double den[3] = {0, 0, 0}, tru[3] = {0, 0, 0}, num[3] = {0, 0, 0}, thm = 0, rep = 0, bd = 0, ag = 0, ag_n = 0;
        for (auto& x : cp.S)
            for (int j = 0; j < 3; j++) {
                double w = x.second * dj[j](x.first);
                den[j] += w, tru[j] += w * std::log2(x.first);
            }
        if (sc.r == 0) {  // context-free: every copy of a selected seed is selected
            for (int j = 0; j < 3; j++) num[j] = tru[j];
            thm = tru[0];
        } else {
            for (auto& x : cp.P.at(sc.r)) {
                long c = (long)x.c, m = (long)x.m;
                double sum = 0, qs = 0;
                for (size_t e = 0; e < x.s.size() && e < q.size(); e++) sum += x.s[e], qs += x.s[e] * q[e];
                double qbar = sum > 0 ? qs / sum : 1, T = Psi(m, c - m, qbar), w = x.positions * d(c);
                for (int j = 0; j < 3; j++) num[j] += x.positions * dj[j](c) * T;
                thm += w * (m == c ? Phi(c - 1, sc.C) : Phi(c - 1, sc.b));
                rep += w, bd += m < c ? w : 0, ag += d(c) * qs, ag_n += d(c) * sum;
            }
        }
        put(row, "H_pred", num[0] / den[0]), put(row, "H_pred_thm", thm / den[0]), put(row, "H_true_pred", tru[0] / den[0]);
        put(row, "H_pred_measured_d", num[1] / den[1]), put(row, "Delta_pred", (tru[0] - num[0]) / den[0]);
        put(row, "H_pred_no_tandem", num[2] / den[2]), put(row, "H_true_pred_no_tandem", tru[2] / den[2]);
        put(row, "beta_emp", sc.r ? bd / rep : 0);
        if (ag_n > 0) put(row, "agree_pred", ag / ag_n);  // mean q over the copies with another context
        double Ltot = 0, beta = 0;
        for (auto& l : cp.E) Ltot += l.second;
        for (auto& l : cp.E) beta += l.second / Ltot * std::min(1.0, 2.0 * sc.r / std::max(l.first - s.k + 1, 1L));
        put(row, "beta_zipf", beta);
        if (tt && !tilt) {  // Zipf law for the aperiodic k-mers, the periodic ones as measured (tandem term)
            Scheme flat = sc;
            flat.d_of_c = nullptr;
            put(row, "H_pred_zipf", zipf_with_tandem(flat, *tt, tandem ? &tandem_calibration(key, s) : nullptr, cp, c_min, N));
        } else {
            put(row, "H_pred_zipf", std::isfinite(alpha) ? predict_zipf(sc, alpha, c_max, cp.E, rho, N) : NaN);
        }
        Scheme flat = sc;  // the earlier form: one Zipf law fitted to distinct seeds, no tandem term
        if (!tilt) flat.d_of_c = nullptr;
        put(row, "H_pred_zipf_distinct", std::isfinite(alpha_d) ? predict_zipf(flat, alpha_d, c_max, cp.E, rho, N) : NaN);
    }
}

// ---------------------------------------------------------------------------
// panel on one genome: for every hash seed, every scheme
// ---------------------------------------------------------------------------
// real genomes: cp (copies tables) and spec (spectrum row) for the predictions, when available
// counted: the genome's k-mer counts at the panel's k, if the caller has them already
static std::vector<Row> panel(const Genome& g, bool real, const CopyTables* cp = nullptr, const Row* spec = nullptr,
                              std::vector<Row>* kernels = nullptr, const Counts* counted = nullptr,
                              const TandemTable* tt = nullptr) {
    Settings s = settings(&g);
    std::vector<std::string> schemes = split(setting(&g, "schemes", PANEL));
    std::vector<std::string> seeds = split(setting(&g, "hash-seeds", real ? "0,1" : "0,1,2"));
    int reads = std::stoi(setting(&g, "reads", real ? "10000" : "0"));
    int read_len = std::stoi(setting(&g, "read-len", real ? "150" : "1000"));
    bool predict_ = has_truth(g) && !real && !opt.count("no-predict");
    if (predict_)  // stop before sketching if a scheme has no calibration
        for (auto& key : schemes)
            if (radius(key, s) > 0) calibration_for(key, s, std::stod(setting(&g, "theta", "0")));

    double t0 = now();
    Counts own_kmers;  // genome k-mer counts, once
    if (!counted) own_kmers = count_space(g, "all", s, 0);
    const Counts& kmers = counted ? *counted : own_kmers;
    std::vector<double> khist = histogram(kmers);
    s.genome = &kmers;
    if (std::find(schemes.begin(), schemes.end(), "oracle") != schemes.end())
        oracle_threshold(kmers, s.d, s.c_star, s.tie);
    if (!counted) log("[%s] counted k=%d in %.1fs", g.name.c_str(), s.k, now() - t0);

    std::vector<uint8_t> in_repeat;  // per genome position: k-mer inside a planted copy
    std::vector<uint32_t> planted;   // per genome position: planted copy number
    if (has_truth(g) && !real) {
        in_repeat.assign(g.offset.back(), 0);
        planted.assign(g.offset.back(), 0);
        for (auto& c : g.copies)
            for (long i = c.start; i + s.k <= c.start + c.length; i++) in_repeat[i] = 1, planted[i] = c.copies;
    }

    std::vector<Row> rows;
    for (auto& hs_str : seeds) {
        uint64_t hs = std::stoull(hs_str);
        log("[%s] hash seed %llu", g.name.c_str(), (unsigned long long)hs);
        std::vector<Row> seed_rows;
        std::map<std::string, std::pair<std::vector<double>, double>> d_of_c;
        for (auto& key : schemes) {
            // real genomes: schemes that do not use the hash seed are measured with the first one only
            if (real && &hs_str != &seeds[0] && (key == "all" || key == "lex-min")) continue;
            double t = now();
            Counts own;  // randstrobes have their own seed space, per hash seed
            if (key == "randstrobe") own = count_space(g, key, s, hs);
            const Counts& sp = key == "randstrobe" ? own : kmers;
            // select, then count every selected seed among the selected (m_S) and in the seed space (m_G)
            std::vector<uint64_t> pos, hash;
            for (size_t r = 0; r < g.seq.size(); r++)
                for (const Seed& x : select_seeds(key, g.seq[r], s, hs)) pos.push_back(g.offset[r] + x.start), hash.push_back(x.hash);
            Counts selected = key == "all" ? Counts() : make_counts(hash);
            std::vector<uint32_t> ms(hash.size()), mg(hash.size());
            std::vector<uint8_t> rep;
            for (size_t i = 0; i < hash.size(); i++) {
                mg[i] = count(sp, hash[i]);
                ms[i] = key == "all" ? mg[i] : count(selected, hash[i]);  // all k-mers: m_S = m_G
                if (!in_repeat.empty()) rep.push_back(in_repeat[pos[i]]);
            }
            Row row;
            put_str(row, "scheme", key);
            Row m = measure(ms, mg, rep.empty() ? nullptr : &rep, sp.sorted.size(), k_eff(key, s), radius(key, s));
            row.insert(row.end(), m.begin(), m.end());
            put_str(row, "scheme_name", scheme_name(key, s));
            put(row, "k", s.k), put(row, "w", s.w), put(row, "hash_seed", (double)hs);
            if (reads > 0) read_stats(row, g, pos, ms, mg, k_eff(key, s), span(key, s), reads, read_len, hs);
            put(row, "seconds", now() - t);
            if (predict_ && (key == "lex-min" || key == "weighted-min")) {
                double overall;
                std::vector<double> rate = selection_function(g, pos, planted, s.k, overall);
                d_of_c[key] = {rate, overall};
            }
            if (cp && key != "randstrobe") {  // real genome: P(selected | genome count c)
                std::vector<double> sel(khist.size(), 0), rate(khist.size(), NaN);
                for (uint32_t c : mg) sel[c]++;
                for (size_t c = 1; c < khist.size(); c++)
                    if (khist[c] > 0) rate[c] = sel[c] / khist[c];
                d_of_c[key] = {rate, (double)hash.size() / kmers.sorted.size()};
            }
            log("    %-40s d=%.4f H=%.4f (%.1fs)", scheme_name(key, s).c_str(), get(row, "d"), get(row, "H"), now() - t);
            seed_rows.push_back(row);
        }
        if (predict_) predictions(seed_rows, g, s, d_of_c);
        if (cp && spec && kernels) real_predictions(seed_rows, s, *cp, *spec, d_of_c, *kernels, tt);
        for (Row& row : seed_rows) {
            put_str(row, "genome", g.name);
            for (auto& kv : g.truth) put(row, kv.first, value(kv.second));
            rows.push_back(row);
        }
    }
    return rows;
}

// ---------------------------------------------------------------------------
// commands
// ---------------------------------------------------------------------------
static int cmd_panel() {
    if (opt.count("calibration")) load_calibration(opt["calibration"]);
    std::vector<Row> rows;
    for (auto& f : files) {
        std::vector<Row> r = panel(load(f), false);
        rows.insert(rows.end(), r.begin(), r.end());
    }
    write_json("panel", {{"rows", rows}});
    return 0;
}

// simulated repeat entropy (log2 genome count of k-mers in planted copies), pooled per alpha
static int cmd_closed_forms() {
    std::map<double, std::vector<double>> per_genome;  // alpha -> mean per genome
    std::map<double, double> num_, den_, fams, cmax;
    std::vector<double> order;
    for (auto& f : files) {
        Genome g = load(f);
        Settings s = settings(&g);
        double a = std::stod(setting(&g, "alpha", "nan"));
        Counts kmers = count_space(g, "all", s, 0);
        std::vector<uint8_t> in_copy(g.seq[0].size(), 0);
        for (auto& c : g.copies)
            for (long i = c.start; i + s.k <= c.start + c.length; i++) in_copy[i] = 1;
        double sum = 0, n = 0;
        for (const Seed& x : all_kmers(g.seq[0], s.k))
            if (in_copy[x.start]) sum += std::log2(count(kmers, x.hash)), n++;
        if (!per_genome.count(a)) order.push_back(a);
        per_genome[a].push_back(sum / n);
        num_[a] += sum, den_[a] += n;
        fams[a] += std::stod(setting(&g, "families", "0"));
        cmax[a] = std::stod(setting(&g, "c_max", "nan"));
    }
    std::vector<Row> rows;
    for (double a : order) {
        ClosedForms cf = closed_forms(a);
        auto& v = per_genome[a];
        double mean = 0, var = 0;
        for (double x : v) mean += x / v.size();
        for (double x : v) var += (x - mean) * (x - mean);
        Row row;
        put(row, "alpha", a), put(row, "H_rep", cf.H_rep), put(row, "H_rep_c2", cf.H_rep_c2);
        put(row, "U_N", cf.U_N), put(row, "U_D", cf.U_D), put(row, "Ehits", cf.Ehits);
        put(row, "H_rep_numeric_cmax1e6", a > 2 ? H_rep_truncated(a, 1000000, 1) : NaN);
        put(row, "c_max_eff", cmax[a]), put(row, "H_rep_c2_truncated", H_rep_truncated(a, (long)cmax[a], 2));
        put(row, "H_rep_sim_c2_pooled", num_[a] / den_[a]);
        put(row, "H_rep_sim_c2_rep_mean", mean), put(row, "H_rep_sim_c2_rep_sd", v.size() > 1 ? std::sqrt(var / (v.size() - 1)) : NaN);
        put(row, "families_total", fams[a]), put(row, "genomes", v.size());
        log("  alpha=%g closed=%.4f truncated=%.4f simulated=%.4f", a, cf.H_rep_c2, get(row, "H_rep_c2_truncated"), num_[a] / den_[a]);
        rows.push_back(row);
    }
    write_json("closed-forms", {{"rows", rows}});
    return 0;
}

// alpha from the genome multiplicities of a FracMinHash sketch; hist = positions by multiplicity
static void alpha_from_hist(Row& row, const std::vector<double>& hist, long c_min) {
    double s = 0, n = 0;
    for (long c = c_min; c < (long)hist.size(); c++) s += hist[c] * std::log2(c), n += hist[c];
    put(row, "H_rep_c2", n ? s / n : NaN);
    put(row, "alpha_hat_H", estimate_alpha(n ? s / n : NaN, true));  // NaN if outside alpha > 2
    put(row, "alpha_hat_mle", fit_zipf_mle(hist, c_min, true, -1));
    put(row, "alpha_hat_mle_trunc", fit_zipf_mle(hist, c_min, true, 0));
}

static int cmd_estimate_alpha() {
    std::vector<Row> rows;
    for (size_t i = 0; i < files.size(); i++) {
        Genome g = load(files[i]);
        Settings s = settings(&g);
        uint64_t hs = std::stoull(split(setting(&g, "hash-seeds", std::to_string(i)))[0]);
        double density = std::stod(setting(&g, "sketch-density", "0.05"));
        Counts kmers = count_space(g, "all", s, 0);
        std::vector<uint64_t> sel;
        for (auto& rec : g.seq)
            for (const Seed& x : fracminhash(rec, s.k, density, hs)) sel.push_back(x.hash);
        Counts sorted = make_counts(sel);
        std::vector<double> hist, hs_hist;
        double H = 0;
        for (uint64_t h : sel) {
            uint32_t mg = count(kmers, h);
            if (hist.size() <= mg) hist.resize(mg + 1, 0);
            hist[mg]++;
            H += std::log2(count(sorted, h));
        }
        Row row;
        put(row, "alpha_true", std::stod(setting(&g, "alpha", "nan"))), put(row, "rep", (double)i);
        put(row, "hash_seed", (double)hs), put(row, "d", (double)sel.size() / kmers.sorted.size());
        put(row, "H_sketch", H / sel.size());
        alpha_from_hist(row, hist, 2);
        put(row, "families", std::stod(setting(&g, "families", "nan")));
        put(row, "c_max_eff", std::stod(setting(&g, "c_max", "nan")));
        put_str(row, "genome", g.name);
        log("  alpha=%g alpha_hat_mle_trunc=%.4f", get(row, "alpha_true"), get(row, "alpha_hat_mle_trunc"));
        rows.push_back(row);
    }
    write_json("estimate-alpha", {{"rows", rows}});
    return 0;
}

// k-mer spectrum and Zipf fits (run_real.spectrum_stats)
static Row spectrum(const std::string& name, int k, const std::vector<double>& hpos, long c_min,
                    Row* spec) {
    std::vector<double> fc(hpos.size(), 0), cs, fs, ps;
    double N = 0, H_all = 0, rep = 0, rep_h = 0, eh = 0, fsum = 0, cmax = 0;
    for (size_t c = 1; c < hpos.size(); c++) {
        fc[c] = std::floor(hpos[c] / c);
        N += hpos[c], H_all += hpos[c] * std::log2(c), eh += hpos[c] * c, fsum += fc[c];
        if (c >= 2) rep += hpos[c], rep_h += hpos[c] * std::log2(c);
        if (fc[c] > 0) cmax = c, cs.push_back(c), fs.push_back(fc[c]), ps.push_back(hpos[c]);
    }
    put_str(*spec, "genome", name), put(*spec, "k", k);
    put(*spec, "c", array(cs)), put(*spec, "f_c", array(fs)), put(*spec, "positions", array(ps));
    double a_trunc = fit_zipf_mle(fc, c_min, false, 0), rho = rep / N;
    double H_rep_zipf = std::isfinite(a_trunc) ? H_rep_truncated(a_trunc, (long)cmax, 2) : NaN;
    // the same law fitted to positions (weight c f_c): its MLE matches the mean of ln c over repeated
    // positions, the quantity H depends on; under one Zipf law both fits agree
    double a_pos = fit_zipf_mle(hpos, c_min, true, 0);
    double H_rep_zipf_pos = std::isfinite(a_pos) ? H_rep_truncated(a_pos, (long)cmax, 2) : NaN;
    double lam = N / std::pow(4.0, k);
    Row r;
    put_str(r, "genome", name), put(r, "k", k), put(r, "positions", N), put(r, "distinct", fsum), put(r, "c_max", cmax);
    put(r, "rho_repeated", rho), put(r, "alpha_mle", fit_zipf_mle(fc, c_min, false, -1)), put(r, "alpha_mle_trunc", a_trunc);
    put(r, "H_all", H_all / N), put(r, "H_rep_measured", rho > 0 ? rep_h / rep : 0), put(r, "H_rep_zipf_trunc", H_rep_zipf);
    put(r, "H_all_zipf_pred", rho * H_rep_zipf), put(r, "lambda_collision", lam), put(r, "H_coll_pred", H_coll(lam));
    put(r, "alpha_pos_trunc", a_pos), put(r, "H_rep_zipf_pos", H_rep_zipf_pos), put(r, "H_all_zipf_pos_pred", rho * H_rep_zipf_pos);
    put(r, "U_N", hpos.size() > 1 ? hpos[1] / N : 0), put(r, "U_D", fc.size() > 1 ? fc[1] / fsum : 0), put(r, "Ehits", eh / N);
    return r;
}

static int cmd_real() {
    std::vector<Row> spec_rows, spectra, panel_rows, ksweep, alpha_rows, copy_rows, extents, kernels, tandem_rows,
        tandem_rate_rows;
    long c_min = std::stol(setting(nullptr, "c-min", "2"));
    double density = std::stod(setting(nullptr, "sketch-density", "0.05"));
    std::string ks = setting(nullptr, "ks", "15,17,19,21,25,31");
    for (auto& f : files) {
        double t0 = now();
        Genome g = load(f);
        g.truth.clear(), g.copies.clear();  // real genome: no predictions or repeat split
        Settings s = settings(&g);
        // 1. spectrum
        Counts kmers = count_space(g, "all", s, 0);
        Row spec;
        Row sp = spectrum(g.name, s.k, histogram(kmers), c_min, &spec);
        spec_rows.push_back(sp), spectra.push_back(spec);
        log("[%s] spectrum: alpha=%.3f rho=%.4f H_all=%.4f", g.name.c_str(), get(sp, "alpha_mle_trunc"),
            get(sp, "rho_repeated"), get(sp, "H_all"));
        // 2. repeat structure from src/copies, for the predictions
        CopyTables cp;
        bool have = !opt.count("no-predict") && load_copies(f, cp);
        if (opt.count("no-predict")) {
        } else if (!have) {
            log("[%s] no %s: no predictions (run src/copies %s %d <radii from src/ue settings> <ks>)",
                g.name.c_str(), cp.path.c_str(), f.c_str(), s.k);
        } else {
            Row cs = copies_summary(g.name, cp, histogram(kmers));
            copy_rows.push_back(cs);
            Row ex;
            std::vector<double> L, nL, Lc, nLc, c, nc;
            for (auto& x : cp.E) L.push_back(x.first), nL.push_back(x.second);
            for (auto& x : cp.C) Lc.push_back(x.first), nLc.push_back(x.second);
            for (auto& x : cp.F) c.push_back(x.first), nc.push_back(x.second);
            put_str(ex, "genome", g.name), put(ex, "k", get(cs, "k"));
            put(ex, "extent", array(L)), put(ex, "extent_positions", array(nL));
            put(ex, "copy_length", array(Lc)), put(ex, "copy_length_copies", array(nLc));
            put(ex, "family_c", array(c)), put(ex, "families", array(nc));
            extents.push_back(ex);
            if (get(cs, "k") != s.k) {
                log("[%s] %s has k=%g: no predictions", g.name.c_str(), cp.path.c_str(), get(cs, "k"));
                have = false;
            } else if (get(cs, "spectrum_mismatch") > 0) {  // another FASTA, or an older copies
                log("[%s] the spectrum in %s differs from ue's in %g positions: no predictions (rerun src/copies)",
                    g.name.c_str(), cp.path.c_str(), get(cs, "spectrum_mismatch"));
                have = false;
            }
            log("[%s] copies: %.0f repeated positions, %.0f copies, %.0f families", g.name.c_str(), get(cs, "repeated"),
                get(cs, "copies"), get(cs, "families"));
        }
        // 2b. tandem term: the periodic k-mers, by period, count and position in their run
        TandemTable tt;
        bool tandem = !opt.count("no-tandem");
        if (tandem) {
            double tt0 = now();
            int dmax = 0;
            for (auto& key : split(setting(&g, "schemes", PANEL))) dmax = std::max(dmax, radius(key, s));
            tt = tandem_table(g, s.k, kmers, dmax);
            double per = 0;
            for (int p = 1; p <= tt.p_max; p++) {
                per += tt.by_p[p];
                Row tr;
                put_str(tr, "genome", g.name), put(tr, "k", s.k), put(tr, "p", p), put(tr, "positions", tt.by_p[p]);
                put(tr, "share", tt.by_p[p] / tt.positions);
                put(tr, "mean_log2_c", tt.by_p[p] > 0 ? tt.logc_by_p[p] / tt.by_p[p] : NaN);
                tandem_rows.push_back(tr);
            }
            put(sp, "tandem_share", per / tt.positions);
            spec_rows.back() = sp;
            log("[%s] tandem: %.0f periodic k-mer positions (%.3f%%) in %.1fs", g.name.c_str(), per,
                100 * per / tt.positions, now() - tt0);
        }
        // 3. scheme panel, with predictions
        std::vector<Row> kern;
        std::vector<Row> r = panel(g, true, have ? &cp : nullptr, &sp, &kern, &kmers, tandem ? &tt : nullptr);
        if (tandem && have)  // run rates of each scheme: interior, and averaged over this genome's runs
            for (auto& key : split(setting(&g, "schemes", PANEL))) {
                if (!tandem_scheme(key)) continue;
                const TandemRates& trs = tandem_calibration(key, s);
                std::vector<double> num(tt.p_max + 1, 0), den(tt.p_max + 1, 0);
                for (auto& x : tt.T) {
                    int p = (x.first >> 16) & 0xff, dl = (x.first >> 8) & 0xff, dr = x.first & 0xff;
                    num[p] += x.second * trs.rel(p, dl, dr), den[p] += x.second;
                }
                for (int p = 1; p <= tt.p_max; p++) {
                    Row tr;
                    put_str(tr, "genome", g.name), put_str(tr, "scheme", key), put(tr, "p", p), put(tr, "r", trs.r);
                    put(tr, "base_rate", trs.base), put(tr, "rel_interior", trs.rel(p, trs.r, trs.r));
                    put(tr, "rel_genome", den[p] > 0 ? num[p] / den[p] : NaN), put(tr, "positions", den[p]);
                    tandem_rate_rows.push_back(tr);
                }
            }
        for (Row& kr : kern) put_str(kr, "genome", g.name), kernels.push_back(kr);
        panel_rows.insert(panel_rows.end(), r.begin(), r.end());
        // 4. alpha from FracMinHash sketches (hash mix64(kmer, 900 + h), as run_real.py)
        for (int h = 0; h < 3; h++) {
            std::vector<double> hist;
            double n = 0;
            for (auto& rec : g.seq)
                for (const Seed& x : fracminhash(rec, s.k, density, 900 + h - 211)) {
                    uint32_t mg = count(kmers, x.hash);
                    if (hist.size() <= mg) hist.resize(mg + 1, 0);
                    hist[mg]++, n++;
                }
            double Hs = 0;
            for (size_t c = 1; c < hist.size(); c++) Hs += hist[c] * std::log2(c);
            Row a;
            put(a, "hash_seed", h), put(a, "d", n / kmers.sorted.size()), put(a, "H_sketch", Hs / n);
            alpha_from_hist(a, hist, c_min);
            put_str(a, "genome", g.name), put(a, "alpha_full_spectrum", get(sp, "alpha_mle_trunc"));
            put(a, "alpha_full_spectrum_pos", get(sp, "alpha_pos_trunc"));
            alpha_rows.push_back(a);
        }
        // 5. k-sweep: H of all k-mers vs k. The spectrum of each k comes from the copies table when it
        // has one (src/copies ks: one pass over the LCP array), the panel's k from step 1, and any
        // other k from counting the genome's k-mers here.
        for (auto& kk : split(ks == "none" ? "" : ks)) {
            Settings sk = s;
            sk.k = std::stoi(kk);
            Row spec_k;
            const char* source = "counted";
            Row st;
            if (sk.k == s.k) {
                st = sp, source = "spectrum above";
            } else if (have && cp.K.count(sk.k)) {
                std::vector<double> hpos(cp.K.at(sk.k).rbegin()->first + 1, 0);
                for (auto& x : cp.K.at(sk.k)) hpos[x.first] = x.second;
                st = spectrum(g.name, sk.k, hpos, c_min, &spec_k), source = "copies";
            } else {
                st = spectrum(g.name, sk.k, histogram(count_space(g, "all", sk, 0)), c_min, &spec_k);
            }
            if (sk.k != s.k) spectra.push_back(spec_k);
            Row row;
            for (const char* key : {"genome", "k", "H_all", "U_N", "Ehits", "alpha_mle_trunc", "lambda_collision", "H_coll_pred"})
                for (auto& kv : st)
                    if (kv.first == key) row.push_back(kv);
            ksweep.push_back(row);
            log("[%s] k=%d: H_all=%.4f (%s)", g.name.c_str(), sk.k, get(st, "H_all"), source);
        }
        log("[%s] done in %.1fs", g.name.c_str(), now() - t0);
    }
    write_json("real", {{"spectrum_summary", spec_rows}, {"spectra", spectra}, {"panel", panel_rows},
                        {"k_sweep", ksweep}, {"alpha_from_sketch", alpha_rows}, {"copies", copy_rows},
                        {"extents", extents}, {"context_kernels", kernels}, {"tandem", tandem_rows},
                        {"tandem_rates", tandem_rate_rows}});
    return 0;
}

// the settings ue uses, and the context radii of the schemes (the radii src/copies needs)
static int cmd_settings() {
    Settings s = settings(nullptr);
    std::set<int> radii;
    for (auto& key : split(setting(nullptr, "schemes", PANEL)))
        if (radius(key, s) > 0) radii.insert(radius(key, s));
    std::string r;
    for (int x : radii) r += (r.empty() ? "" : ",") + std::to_string(x);
    printf("k\t%d\nw\t%d\ndensity\t%s\nmodmin_w\t%d\nradii\t%s\n", s.k, s.w, num(s.d).c_str(), s.modmin_w,
           r.empty() ? "0" : r.c_str());
    return 0;
}

static int cmd_sketch() {
    std::string key = setting(nullptr, "scheme", "random-min");
    uint64_t hs = std::stoull(split(setting(nullptr, "hash-seeds", "0"))[0]);
    printf("genome\trecord\tstart\tend\thash\n");
    for (auto& f : files) {
        Genome g = load(f);
        Settings s = settings(&g);
        Counts kmers;
        if (key == "weighted-min" || key == "oracle") {
            kmers = count_space(g, "all", s, 0);
            s.genome = &kmers;
            oracle_threshold(kmers, s.d, s.c_star, s.tie);
        }
        log("[%s] %s", g.name.c_str(), scheme_name(key, s).c_str());
        for (size_t r = 0; r < g.seq.size(); r++)
            for (const Seed& x : select_seeds(key, g.seq[r], s, hs))
                printf("%s\t%s\t%u\t%u\t%llu\n", g.name.c_str(), g.names[r].c_str(), x.start, x.end, (unsigned long long)x.hash);
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2 || !strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) {
        fputs(USAGE, stdout);
        return argc < 2;
    }
    try {
        for (int i = 2; i < argc; i++) {
            std::string a = argv[i];
            if (a.rfind("--", 0) != 0) {
                files.push_back(a);
                continue;
            }
            std::string key = a.substr(2);
            if (FLAGS.count(key)) opt[key] = "1";
            else if (VALUES.count(key) && i + 1 < argc) opt[key] = argv[++i];
            else throw std::runtime_error("unknown option or missing value: " + a);
        }
        std::string c = argv[1];
        if (opt.count("help")) return fputs(USAGE, stdout), 0;
        if (c == "panel") return cmd_panel();
        if (c == "closed-forms") return cmd_closed_forms();
        if (c == "estimate-alpha") return cmd_estimate_alpha();
        if (c == "real") return cmd_real();
        if (c == "sketch") return cmd_sketch();
        if (c == "settings") return cmd_settings();
        throw std::runtime_error("unknown command " + c);
    } catch (const std::exception& e) {
        fprintf(stderr, "ue: %s\n", e.what());
        return 2;
    }
}
