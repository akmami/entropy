// Theorem 1 predictions, Zipf closed forms, alpha estimators, and the empirical
// calibration of C and b. A port of ue/theory.py and ue/calibrate.py.
#pragma once

#include <algorithm>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "sketches.h"

// ---- Theorem 1 --------------------------------------------------------------
typedef std::vector<std::pair<long, long>> Families;  // (copies c, length L), one per family

// How a scheme enters Theorem 1 (ue.theory.SchemeParams)
struct Scheme {
    int k = 21;                             // k_S
    int r = 0;                              // context radius
    double C = 1, b = 1;                    // interior co-selection, boundary agreement
    double coll_agree = 1;                  // co-selection of a random background collision
    std::function<double(double)> b_of_L;   // optional: b by copy length
    std::function<double(long)> d_of_c;     // optional: selection function (lex/weighted minimizers)
};

struct Prediction { double H, H_genome, Delta_order, Delta_cons, Delta_ctx; };

// m = n - k + 1 positions; n = genome length, which sets the collision rate n / 4^k.
// outer_survival = shared_mutations = true, edges = false is Eq. (master) of the paper; edges = true (the
// default, used in the experiments) also counts the k-mers that straddle copy edges. All false is an earlier form.
Prediction predict(const Scheme& s, const Families& fam, double m, double n, double theta, bool outer_survival = true, bool edges = true, bool shared_mutations = true);

// Theorem 1 under a Zipf law fitted to a real genome (theta = 0): copy numbers on [2, c_max] with
// position weight c^(1-alpha), copy lengths (bases) with position weights `lengths`, a fraction rho
// of positions repeated.
double predict_zipf(const Scheme& s, double alpha, long c_max, const std::map<long, double>& lengths, double rho, double n);

double H_coll(double lam);            // E[log2(1 + Poisson(lam))]
double Phi(long n_sib, double q);     // E[log2(1 + Bin(n_sib, q))]
double Psi(long a, long n, double q);  // E[log2(a + Bin(n, q))]
double crude_b(double d);
double crude_C(double theta, int w, double d);

// ---- Zipf closed forms and alpha ---------------------------------------------
double zeta(double s, double q = 1);  // Hurwitz zeta
double H_rep_closed(double alpha);    // c >= 1
double H_rep_c2(double alpha);        // c >= 2
double H_rep_truncated(double alpha, long c_max, long c_min = 2);
struct ClosedForms { double H_rep, H_rep_c2, U_N, U_D, Ehits; };
ClosedForms closed_forms(double alpha);
double estimate_alpha(double H_rep, bool c_ge_2);  // NaN if out of range
// discrete power-law MLE from hist[c]; c_max: -1 untruncated, 0 largest observed
double fit_zipf_mle(const std::vector<double>& hist, long c_min, bool position_weighted, long c_max);

// ---- calibration --------------------------------------------------------------
typedef std::function<std::vector<Seed>(const std::string&)> Sketch;
// the scheme as calibration runs it: hash seed 0; the weighted minimizer counts the sequence itself
Sketch calibration_sketch(const std::string& key, const Settings& s);

struct Calibration {
    double p, C, b, d;
    std::map<int, double> b_profile;  // boundary agreement by distance to the copy edge
    long n_interior, n_boundary;
};

// select: the scheme; space: its seed space (all k-mers, or unthinned randstrobes);
// r: context radius; span: bases one seed position needs
Calibration calibrate(Sketch select, Sketch space, int r, int span, int k, double theta, int L, int reps, uint64_t seed = 0);
std::function<double(double)> b_of_length(Sketch select, Sketch space, int r, int span, int k, double theta, std::vector<int> lengths, int reps, uint64_t seed = 0);
double collision_agreement(Sketch select, int r, int k, int reps, uint64_t seed = 0);
// Agreement between two copies of a seed whose contexts agree on exactly e bases on either side
// and are independent beyond: q[e] = P[the copy is selected | the seed is selected], e = 0..r-1,
// each from `hits` selected seeds. (q[0] is the agreement of a chance collision.)
std::vector<double> context_kernel(Sketch select, int r, int k, long hits = 10000, uint64_t seed = 0);

// ---- tandem term ----------------------------------------------------------------
// Selection rate of a periodic k-mer (minimal period p <= p_max) inside a tandem run, relative to the
// scheme's rate on aperiodic sequence (base), averaged over random units and hash seeds. The run is
// the stretch of consecutive periodic k-mers; dl, dr are the k-mer's distances to its two ends, capped
// at r (r = 0 for context-free schemes: one rate per period). rel(p, r, r) is the interior rate,
// 1/(p d) for leftmost-argmin schemes with p up to the window (Lemma: ties in tandem repeats).
typedef std::function<std::vector<Seed>(const std::string&, uint64_t)> SketchHS;  // sequence, hash seed
struct TandemRates {
    int r = 0, p_max = 0;
    double base = 0;                 // selection rate on aperiodic sequence
    std::vector<double> rate, n;     // [(p * (r + 1) + dl) * (r + 1) + dr]
    double rel(int p, int dl, int dr) const {
        if (p < 1 || p > p_max) return 1;
        size_t i = ((size_t)p * (r + 1) + std::min(dl, r)) * (r + 1) + std::min(dr, r);
        return base > 0 ? rate[i] / base : 1;
    }
};
TandemRates tandem_rates(SketchHS select, int r, int k, int p_max, int reps_per_length = 150, uint64_t seed = 0);
