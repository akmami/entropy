// All seed selection schemes; each mirrors the class of the same name in ue/schemes.py.
#include "sketches.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>

static const uint64_t MAX64 = ~0ULL;

uint64_t mix64(uint64_t x, uint64_t seed) {
    uint64_t z = x + 0x9E3779B97F4A7C15ULL * (seed + 1);
    z ^= z >> 30;
    z *= 0xBF58476D1CE4E5B9ULL;
    z ^= z >> 27;
    z *= 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static uint64_t combine(uint64_t a, uint64_t b, uint64_t seed) {
    return mix64((mix64(a, seed) * 0x100000001B3ULL) ^ b, seed + 7);
}

// hash < f * 2^64, as min(int(f * 2**64), 2**64 - 1) in Python
static uint64_t below(double f) {
    double x = f * 18446744073709551616.0;
    if (x <= 0) return 0;
    if (x >= 18446744073709551616.0) return MAX64;
    return (uint64_t)x;
}

static uint8_t code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

// Split seq at non-ACGT bases; call fn(offset, codes) for each piece with >= min_len bases.
template <class F>
static void for_each_piece(const std::string& seq, size_t min_len, F fn) {
    std::vector<uint8_t> c;
    size_t i = 0;
    while (i < seq.size()) {
        while (i < seq.size() && code(seq[i]) == 4) i++;
        size_t start = i;
        c.clear();
        while (i < seq.size() && code(seq[i]) != 4) c.push_back(code(seq[i++]));
        if (!c.empty() && c.size() >= min_len) fn(start, c);
    }
}

// Packed k-mers of a piece (2 bits per base, forward strand).
static void kmers(const std::vector<uint8_t>& c, int k, std::vector<uint64_t>& v) {
    v.clear();
    uint64_t mask = k == 32 ? MAX64 : (1ULL << (2 * k)) - 1, fwd = 0;
    for (size_t i = 0; i < c.size(); i++) {
        fwd = ((fwd << 2) | c[i]) & mask;
        if (i + 1 >= (size_t)k) v.push_back(fwd);
    }
}

// For each window of W consecutive values, the index of its leftmost minimum.
static std::vector<uint32_t> window_argmin(const std::vector<uint64_t>& x, size_t W) {
    std::vector<uint32_t> out, q(x.size());  // q: monotone queue of indices
    size_t head = 0, tail = 0;
    for (size_t j = 0; j < x.size(); j++) {
        while (tail > head && x[q[tail - 1]] > x[j]) tail--;
        q[tail++] = (uint32_t)j;
        if (j + 1 >= W) {
            while (q[head] + W <= j) head++;
            out.push_back(q[head]);
        }
    }
    return out;
}

static Seed kmer_seed(uint64_t kmer, size_t pos, int k) {
    return {mix64(kmer, 0), (uint32_t)pos, (uint32_t)(pos + k)};
}

// ---------------------------------------------------------------------------
// context-free schemes: a k-mer is selected by its own content
// ---------------------------------------------------------------------------
std::vector<Seed> all_kmers(const std::string& seq, int k) {
    std::vector<Seed> out;
    std::vector<uint64_t> v;
    for_each_piece(seq, k, [&](size_t off, const std::vector<uint8_t>& c) {
        kmers(c, k, v);
        for (size_t i = 0; i < v.size(); i++) out.push_back(kmer_seed(v[i], off + i, k));
    });
    return out;
}

std::vector<Seed> fracminhash(const std::string& seq, int k, double f, uint64_t seed) {
    std::vector<Seed> out;
    std::vector<uint64_t> v;
    for_each_piece(seq, k, [&](size_t off, const std::vector<uint8_t>& c) {
        kmers(c, k, v);
        for (size_t i = 0; i < v.size(); i++)
            if (mix64(v[i], seed + 211) < below(f)) out.push_back(kmer_seed(v[i], off + i, k));
    });
    return out;
}

std::vector<Seed> modmer(const std::string& seq, int k, uint64_t m, uint64_t seed) {
    std::vector<Seed> out;
    std::vector<uint64_t> v;
    for_each_piece(seq, k, [&](size_t off, const std::vector<uint8_t>& c) {
        kmers(c, k, v);
        for (size_t i = 0; i < v.size(); i++)
            if (mix64(v[i], seed + 307) % m == 0) out.push_back(kmer_seed(v[i], off + i, k));
    });
    return out;
}

// A k-mer is a syncmer by the offset (0..k-s) of its smallest s-mer.
static std::vector<Seed> syncmer(const std::string& seq, int k, int s, uint64_t seed, bool open) {
    std::vector<Seed> out;
    std::vector<uint64_t> v, sh;
    for_each_piece(seq, k, [&](size_t off, const std::vector<uint8_t>& c) {
        kmers(c, k, v);
        kmers(c, s, sh);
        for (auto& x : sh) x = mix64(x, seed + 101);
        std::vector<uint32_t> arg = window_argmin(sh, k - s + 1);  // one window per k-mer
        for (size_t i = 0; i < v.size(); i++) {
            int o = (int)(arg[i] - i);
            bool sel = open ? o == (k - s) / 2 : (o == 0 || o == k - s);
            if (sel) out.push_back(kmer_seed(v[i], off + i, k));
        }
    });
    return out;
}

std::vector<Seed> open_syncmer(const std::string& seq, int k, int s, uint64_t seed) {
    return syncmer(seq, k, s, seed, true);
}

std::vector<Seed> closed_syncmer(const std::string& seq, int k, int s, uint64_t seed) {
    return syncmer(seq, k, s, seed, false);
}

Counts make_counts(std::vector<uint64_t> hashes) {
    Counts c;
    std::sort(hashes.begin(), hashes.end());
    c.sorted = std::move(hashes);
    int bits = 1;
    while (bits < 24 && (4ULL << bits) < c.sorted.size()) bits++;  // about 4 entries per index slot
    c.shift = 64 - bits;
    c.index.assign((1ULL << bits) + 1, c.sorted.size());
    for (size_t i = c.sorted.size(); i-- > 0;) c.index[c.sorted[i] >> c.shift] = i;
    for (size_t b = c.index.size() - 1; b-- > 0;) c.index[b] = std::min(c.index[b], c.index[b + 1]);
    return c;
}

uint32_t count(const Counts& c, uint64_t hash) {
    if (c.sorted.empty()) return 0;
    auto first = c.sorted.begin() + c.index[hash >> c.shift], last = c.sorted.begin() + c.index[(hash >> c.shift) + 1];
    auto r = std::equal_range(first, last, hash);
    return (uint32_t)(r.second - r.first);
}

void oracle_threshold(const Counts& counts, double d, uint32_t& c_star, double& tie_fraction) {
    const std::vector<uint64_t>& genome = counts.sorted;
    std::vector<double> hist;  // positions per genome count
    for (size_t i = 0; i < genome.size();) {
        size_t j = i;
        while (j < genome.size() && genome[j] == genome[i]) j++;
        if (hist.size() <= j - i) hist.resize(j - i + 1, 0);
        hist[j - i] += (double)(j - i);
        i = j;
    }
    double target = d * (double)genome.size(), cum = 0;
    c_star = (uint32_t)hist.size();
    tie_fraction = 0;
    for (size_t c = 0; c < hist.size(); c++) {
        if (cum + hist[c] >= target) {
            c_star = (uint32_t)c;
            if (hist[c] > 0) tie_fraction = (target - cum) / hist[c];
            return;
        }
        cum += hist[c];
    }
}

std::vector<Seed> oracle(const std::string& seq, int k, const Counts& genome, uint32_t c_star, double tie_fraction, uint64_t seed) {
    std::vector<Seed> out;
    std::vector<uint64_t> v;
    for_each_piece(seq, k, [&](size_t off, const std::vector<uint8_t>& c) {
        kmers(c, k, v);
        for (size_t i = 0; i < v.size(); i++) {
            uint32_t gc = count(genome, mix64(v[i], 0));
            bool sel = gc > 0 && gc < c_star;
            if (gc == c_star) sel = mix64(v[i], seed + 601) < below(tie_fraction);
            if (sel) out.push_back(kmer_seed(v[i], off + i, k));
        }
    });
    return out;
}

std::vector<Seed> randstrobe(const std::string& seq, int ks, int w_min, int w_max, uint64_t seed) {
    std::vector<Seed> out;
    std::vector<uint64_t> s, h;
    for_each_piece(seq, ks + w_max, [&](size_t off, const std::vector<uint8_t>& c) {
        kmers(c, ks, s);
        h.resize(s.size());
        for (size_t a = 0; a < s.size(); a++) h[a] = mix64(s[a], seed + 503);
        for (size_t a = 0; a + w_max < s.size(); a++) {
            // second strobe: the j in [a + w_min, a + w_max] minimising (h[a] + h[j]) mod 2^64
            uint64_t best = MAX64;
            size_t j = a;
            for (size_t b = a + w_min; b <= a + w_max; b++) {
                uint64_t sc = h[b] == MAX64 ? MAX64 : h[a] + h[b];
                if (sc < best) best = sc, j = b;
            }
            if (h[a] == MAX64 || best == MAX64) continue;
            out.push_back({combine(s[a], s[j], 0), (uint32_t)(off + a), (uint32_t)(off + j + ks)});
        }
    });
    return out;
}

// ---------------------------------------------------------------------------
// context-dependent schemes: the smallest k-mer of every window of w k-mers
// ---------------------------------------------------------------------------
template <class Key>
static std::vector<Seed> minimizer(const std::string& seq, int k, int w, Key key) {
    std::vector<Seed> out;
    std::vector<uint64_t> v, key_of;
    std::vector<uint8_t> mark;
    for_each_piece(seq, w + k - 1, [&](size_t off, const std::vector<uint8_t>& c) {
        kmers(c, k, v);
        key_of.resize(v.size());
        for (size_t i = 0; i < v.size(); i++) key_of[i] = key(v[i]);
        mark.assign(v.size(), 0);
        for (uint32_t i : window_argmin(key_of, w)) mark[i] = 1;
        for (size_t i = 0; i < v.size(); i++)
            if (mark[i]) out.push_back(kmer_seed(v[i], off + i, k));
    });
    return out;
}

std::vector<Seed> random_minimizer(const std::string& seq, int k, int w, uint64_t seed) {
    return minimizer(seq, k, w, [&](uint64_t x) { return mix64(x, seed); });
}

std::vector<Seed> lex_minimizer(const std::string& seq, int k, int w) {
    return minimizer(seq, k, w, [](uint64_t x) { return x; });
}

std::vector<Seed> weighted_minimizer(const std::string& seq, int k, int w, uint64_t seed, const Counts& genome, uint32_t max_count) {
    return minimizer(seq, k, w, [&](uint64_t x) {
        bool frequent = count(genome, mix64(x, 0)) > max_count;
        return (mix64(x, seed) >> 1) | (frequent ? 1ULL << 63 : 0);
    });
}

// Mod-minimizer (Groot Koerkamp & Pibiri 2024): t = 4 + ((k-4) mod w); the smallest
// t-mer x of the window (w+k-t t-mers) selects the k-mer at offset x mod w.
std::vector<Seed> mod_minimizer(const std::string& seq, int k, int w, uint64_t seed) {
    int t = 4 + (k - 4) % w;
    std::vector<Seed> out;
    std::vector<uint64_t> v, th;
    std::vector<uint8_t> mark;
    for_each_piece(seq, w + k - 1, [&](size_t off, const std::vector<uint8_t>& c) {
        kmers(c, k, v);
        kmers(c, t, th);
        for (auto& x : th) x = mix64(x, seed + 701);
        std::vector<uint32_t> arg = window_argmin(th, w + k - t);
        mark.assign(v.size(), 0);
        for (size_t i = 0; i < arg.size(); i++) mark[i + (arg[i] - i) % w] = 1;
        for (size_t i = 0; i < v.size(); i++)
            if (mark[i]) out.push_back(kmer_seed(v[i], off + i, k));
    });
    return out;
}

std::vector<PeriodicKmer> periodic_kmers(const std::string& seq, int k, int p_max) {
    std::vector<PeriodicKmer> out;
    std::vector<uint64_t> v;
    p_max = std::min(p_max, k / 2);
    for_each_piece(seq, k, [&](size_t off, const std::vector<uint8_t>& c) {
        kmers(c, k, v);
        for (size_t i = 0; i < v.size(); i++)
            for (int p = 1; p <= p_max; p++) {  // packed: first base most significant
                uint64_t tail = (1ULL << (2 * (k - p))) - 1;
                if ((v[i] >> (2 * p)) == (v[i] & tail)) {
                    out.push_back({mix64(v[i], 0), (uint32_t)(off + i), (uint8_t)p});
                    break;
                }
            }
    });
    return out;
}

std::vector<Seed> thin(const std::vector<Seed>& seeds, double f, uint64_t seed) {
    if (f >= 1) return seeds;
    std::vector<Seed> out;
    for (const Seed& s : seeds)
        if (mix64(s.hash, seed + 401) < below(f)) out.push_back(s);
    return out;
}

// ---------------------------------------------------------------------------
// the panel: every scheme at matched density d (ue.experiment.make_scheme)
// ---------------------------------------------------------------------------
static int open_s(int k, double d) { return std::max(std::min(8, k - 2), k + 1 - (int)std::floor(1.0 / d)); }
static int closed_s(int k, double d) { return std::max(std::min(8, k - 2), (int)std::ceil(k + 1 - 2.0 / d)); }
static int strobe(int k) { return (k + 1) / 2; }  // strobe length: k_eff = 2 * strobe >= k

std::vector<Seed> select_seeds(const std::string& key, const std::string& seq, const Settings& s, uint64_t hs) {
    int k = s.k;
    double d = s.d;
    if (key == "all") return all_kmers(seq, k);
    if (key == "random-min") return random_minimizer(seq, k, s.w, hs);
    if (key == "lex-min") return lex_minimizer(seq, k, s.w);
    if (key == "weighted-min") return weighted_minimizer(seq, k, s.w, hs, *s.genome, 1);
    if (key == "mod-min") return mod_minimizer(seq, k, s.modmin_w, hs);
    if (key == "open-sync") {
        int sm = open_s(k, d);
        return thin(open_syncmer(seq, k, sm, hs), std::min(1.0, d / (1.0 / (k - sm + 1))), hs);
    }
    if (key == "closed-sync") {
        int sm = closed_s(k, d);
        return thin(closed_syncmer(seq, k, sm, hs), std::min(1.0, d / (2.0 / (k - sm + 1))), hs);
    }
    if (key == "fracminhash") return fracminhash(seq, k, d, hs);
    if (key == "modmer") return modmer(seq, k, (uint64_t)std::max(1.0, std::nearbyint(1.0 / d)), hs);
    if (key == "randstrobe") return thin(randstrobe(seq, strobe(k), strobe(k) + 1, strobe(k) + 40, hs), d, hs);
    if (key == "oracle") return oracle(seq, k, *s.genome, s.c_star, s.tie, hs);
    throw std::runtime_error("unknown scheme " + key);
}

std::vector<Seed> seed_space(const std::string& key, const std::string& seq, const Settings& s, uint64_t hs) {
    if (key == "randstrobe") return randstrobe(seq, strobe(s.k), strobe(s.k) + 1, strobe(s.k) + 40, hs);
    return all_kmers(seq, s.k);
}

int radius(const std::string& key, const Settings& s) {
    if (key == "mod-min") return s.modmin_w - 1;
    return key == "random-min" || key == "lex-min" || key == "weighted-min" ? s.w - 1 : 0;
}
int k_eff(const std::string& key, const Settings& s) { return key == "randstrobe" ? 2 * strobe(s.k) : s.k; }
int span(const std::string& key, const Settings& s) {
    return key == "randstrobe" ? 2 * strobe(s.k) + 40 : s.k;  // bases one seed position needs
}

std::string scheme_name(const std::string& key, const Settings& s) {  // as in misc/ue/schemes.py
    char b[96];
    int k = s.k;
    double d = s.d;
    if (key == "all") return "all";
    if (key == "random-min" || key == "lex-min" || key == "weighted-min")
        snprintf(b, sizeof b, "%s-minimizer(w=%d)", key.substr(0, key.size() - 4).c_str(), s.w);
    else if (key == "mod-min") snprintf(b, sizeof b, "mod-minimizer(w=%d)", s.modmin_w);
    else if (key == "open-sync") snprintf(b, sizeof b, "open-syncmer(s=%d)", open_s(k, d));
    else if (key == "closed-sync") snprintf(b, sizeof b, "closed-syncmer(s=%d)", closed_s(k, d));
    else if (key == "fracminhash") snprintf(b, sizeof b, "fracminhash(f=%.4g)", d);
    else if (key == "modmer") snprintf(b, sizeof b, "modmer(m=%.0f)", std::max(1.0, std::nearbyint(1.0 / d)));
    else if (key == "randstrobe") snprintf(b, sizeof b, "randstrobe(2,%d,%d,%d)", strobe(k), strobe(k) + 1, strobe(k) + 40);
    else snprintf(b, sizeof b, "oracle(d=%.4g)", d);
    return b;
}

// mod-minimizer window whose density on a random string is closest to d
// (the search of misc/ue/experiment.py's _modmin_w_for_density, on a C++ random
// string)
int modmin_window(int k, double d) {
    std::mt19937_64 rng(12345);
    std::string seq(200000, 'A');
    for (auto& c : seq) c = "ACGT"[rng() >> 62];
    int best = 2;
    double best_err = 1e9;
    for (int w = 2; w < 8 * (int)std::nearbyint(1 / d) + 8; w++) {
        double dm = (double)mod_minimizer(seq, k, w, 0).size() / (seq.size() - k + 1);
        if (std::fabs(dm - d) < best_err) best = w, best_err = std::fabs(dm - d);
        if (dm < d) break;
    }
    return best;
}
