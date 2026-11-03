// All seed selection schemes. Each one takes a sequence (one FASTA record)
// and returns its seeds in the same struct. Bases other than ACGT/acgt split
// the sequence: no seed and no window crosses them.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

struct Seed {
    uint64_t hash;   // seed identity: equal hash <=> equal seed
    uint32_t start;  // first base, 0-based
    uint32_t end;    // one past the last base
};

// Genome-wide seed counts: the sorted hashes of every seed of a genome, plus an
// index on their top bits so a lookup searches only a few entries.
struct Counts {
    std::vector<uint64_t> sorted;
    std::vector<uint64_t> index;  // index[b] = first entry whose top bits are >= b
    int shift = 64;
};
Counts make_counts(std::vector<uint64_t> hashes);
uint32_t count(const Counts& c, uint64_t hash);

// k-mer schemes: hash = mix64(packed k-mer, 0), a bijection. Every scheme reads the forward
// strand only: a k-mer and its reverse complement are different seeds.
// `seed` is the hash seed; schemes use it exactly as ue/schemes.py does.
std::vector<Seed> all_kmers(const std::string& seq, int k);
std::vector<Seed> random_minimizer(const std::string& seq, int k, int w, uint64_t seed);
std::vector<Seed> lex_minimizer(const std::string& seq, int k, int w);
// k-mers with genome count > max_count are ordered after all others (Winnowmap)
std::vector<Seed> weighted_minimizer(const std::string& seq, int k, int w, uint64_t seed, const Counts& genome, uint32_t max_count);
std::vector<Seed> mod_minimizer(const std::string& seq, int k, int w, uint64_t seed);
std::vector<Seed> open_syncmer(const std::string& seq, int k, int s, uint64_t seed);
std::vector<Seed> closed_syncmer(const std::string& seq, int k, int s, uint64_t seed);
std::vector<Seed> fracminhash(const std::string& seq, int k, double f, uint64_t seed);
std::vector<Seed> modmer(const std::string& seq, int k, uint64_t m, uint64_t seed);
// positions of lowest genome count: count < c_star, plus a hash-chosen fraction of count == c_star
std::vector<Seed> oracle(const std::string& seq, int k, const Counts& genome, uint32_t c_star, double tie_fraction, uint64_t seed);
// order-2 randstrobes; hash = combine(strobe1, strobe2, 0)
std::vector<Seed> randstrobe(const std::string& seq, int ks, int w_min, int w_max, uint64_t seed);

// keep a hash-chosen fraction f of the seeds (FracMinHash thinning, to match densities)
std::vector<Seed> thin(const std::vector<Seed>& seeds, double f, uint64_t seed);

// k-mers that are periodic: minimal period p <= p_max (p_max <= k/2), i.e. G[i+j] = G[i+j+p] for
// 0 <= j < k-p. Inside a tandem run every k-mer is periodic with the run's period. hash as all_kmers.
struct PeriodicKmer { uint64_t hash; uint32_t start; uint8_t p; };
std::vector<PeriodicKmer> periodic_kmers(const std::string& seq, int k, int p_max);

// c_star and tie_fraction of the oracle at density d
void oracle_threshold(const Counts& genome, double d, uint32_t& c_star, double& tie_fraction);

uint64_t mix64(uint64_t x, uint64_t seed);  // splitmix64, as ue.seq.mix64

// ---- the panel: every scheme at a matched density d (ue, calibrate) ----------------------------
// keys: all, random-min, lex-min, weighted-min, mod-min, open-sync, closed-sync, fracminhash, modmer,
// randstrobe, oracle
struct Settings {
    int k = 21, w = 11, modmin_w = 0;     // k, window, mod-minimizer window (searched for d)
    double d = 0;                         // target density, 2 / (w + 1)
    const Counts* genome = nullptr;       // counts of all k-mers (weighted minimizer, oracle)
    uint32_t c_star = 0;                  // oracle threshold
    double tie = 0;
};
// the scheme's seeds at density s.d; syncmers and randstrobes are thinned to it
std::vector<Seed> select_seeds(const std::string& key, const std::string& seq, const Settings& s, uint64_t hash_seed);
// the seed space its genome multiplicities m_G are counted in (all k-mers; unthinned randstrobes)
std::vector<Seed> seed_space(const std::string& key, const std::string& seq, const Settings& s, uint64_t hash_seed);
int radius(const std::string& key, const Settings& s);  // context radius r
int k_eff(const std::string& key, const Settings& s);   // effective seed length
int span(const std::string& key, const Settings& s);    // bases one seed position needs
std::string scheme_name(const std::string& key, const Settings& s);
int modmin_window(int k, double d);  // mod-minimizer window whose density on a random string is closest to d
