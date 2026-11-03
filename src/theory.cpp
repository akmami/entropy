// Port of ue/theory.py and ue/calibrate.py; names follow the Python code.
#include "theory.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <tuple>

static const double LN2 = std::log(2.0);
static const double NaN = std::nan("");
static const double Inf = INFINITY;

// ---------------------------------------------------------------------------
// building blocks
// ---------------------------------------------------------------------------
double H_coll(double lam) {
    if (lam <= 0) return 0;
    if (lam > 1e4) return std::log2(1 + lam) - lam / (2 * (1 + lam) * (1 + lam) * LN2);
    long jmax = (long)(lam + 12 * std::sqrt(lam) + 30);
    double lp = -lam, s = 0;  // log Poisson pmf, by recursion
    for (long j = 0; j <= jmax; j++) {
        if (j > 0) lp += std::log(lam) - std::log((double)j);
        s += std::exp(lp) * std::log2(1.0 + j);
    }
    return s;
}

double Phi(long n, double q) {
    if (n <= 0 || q <= 0) return 0;
    if (q >= 1) return std::log2(1.0 + n);
    static std::map<std::pair<long, double>, double> cache;
    auto it = cache.find({n, q});
    if (it != cache.end()) return it->second;
    double mu = n * q, var = n * q * (1 - q), s = 0;
    if (n > 3000 && mu > 50) {
        s = std::log2(1 + mu) - var / (2 * (1 + mu) * (1 + mu) * LN2);
    } else {
        double lp = n * std::log1p(-q);  // log binomial pmf, by recursion
        for (long j = 0; j <= n; j++) {
            if (j > 0) lp += std::log((double)(n - j + 1) / j) + std::log(q) - std::log1p(-q);
            s += std::exp(lp) * std::log2(1.0 + j);
        }
    }
    return cache[{n, q}] = s;
}

// E[log2(a + Bin(n, q))]: a copies selected for sure, n more with probability q each
double Psi(long a, long n, double q) {
    if (n <= 0 || q <= 0) return std::log2((double)a);
    if (q >= 1) return std::log2((double)(a + n));
    static std::map<std::tuple<long, long, double>, double> cache;
    auto key = std::make_tuple(a, n, q);
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    double mu = n * q, var = n * q * (1 - q), s = 0;
    if (n > 3000 && mu > 50) {  // delta method
        s = std::log2(a + mu) - var / (2 * (a + mu) * (a + mu) * LN2);
    } else {  // the binomial pmf within 12 sd of the mean, by recursion from the lowest term
        double sd = std::sqrt(var);
        long lo = std::max(0L, (long)std::floor(mu - 12 * sd - 5)), hi = std::min(n, (long)std::ceil(mu + 12 * sd + 5));
        double lp = std::lgamma(n + 1.0) - std::lgamma(lo + 1.0) - std::lgamma(n - lo + 1.0) + lo * std::log(q) + (n - lo) * std::log1p(-q);
        double p = std::exp(lp), ratio = q / (1 - q);
        for (long j = lo; j <= hi; j++) {
            s += p * std::log2((double)(a + j));
            p *= (double)(n - j) / (j + 1) * ratio;
        }
    }
    return cache[key] = s;
}

double crude_b(double d) { return (1 + d) / 2; }

double crude_C(double theta, int w, double d) {
    double p_ctx = std::pow(1 - theta, w - 1);
    return p_ctx + (1 - p_ctx) * crude_b(d);
}

static double binom_pmf(int j, int n, double p) {
    if (p <= 0) return j == 0;
    return std::exp(std::lgamma(n + 1.0) - std::lgamma(j + 1.0) - std::lgamma(n - j + 1.0) + j * std::log(p) + (n - j) * std::log1p(-p));
}

// ---------------------------------------------------------------------------
// Theorem 1 (ue.theory._expected_H)
// ---------------------------------------------------------------------------
static double expected_H(const Scheme& s, const Families& fam, double m, double n, double theta,
                         bool outer_survival, bool edges, bool shared_mutations, double C,
                         const std::function<double(double)>& b_of_L, bool tilt, double coll_agree) {
    int k = s.k;
    double p = std::pow(1 - theta, k), outer = outer_survival ? p : 1;
    bool use_tilt = tilt && s.d_of_c;
    double wsum = 0, hsum = 0, used = 0;
    for (auto& f : fam) {
        long c = f.first, L = f.second;
        double dc = use_tilt ? s.d_of_c(c) : 1;
        double n_int = (double)c * std::max(L - k + 1, 0L);
        double beta = std::min(1.0, 2.0 * s.r / std::max(L - k + 1, 1L));
        double b = b_of_L(L), T = 0;
        if (outer_survival && theta > 0 && shared_mutations) {
            // own seed carries j mutations; a sibling carries the identical seed with prob (1-θ)^(k-j) (θ/3)^j
            for (int j = 0; j <= std::min(k, 4); j++) {
                double qj = std::pow(1 - theta, k - j) * std::pow(theta / 3, j);
                T += binom_pmf(j, k, theta) * ((1 - beta) * Phi(c - 1, qj * C) + beta * Phi(c - 1, qj * b));
            }
        } else {
            T = outer * ((1 - beta) * Phi(c - 1, p * C) + beta * Phi(c - 1, p * b));
        }
        wsum += n_int * dc;
        hsum += n_int * dc * T;
        used += n_int;
        if (edges && c > 1) {  // seeds spanning a copy edge with j random flank bases
            double n_edge = 2.0 * c * (k - 1), E = 0, ag = s.r > 0 ? b : 1;
            for (int j = 1; j < k; j++) {
                double pk = std::pow(1 - theta, k - j);
                E += (outer_survival ? pk : 1) * Phi(c - 1, pk * std::ldexp(1.0, -2 * j) * ag);
            }
            E /= k - 1;
            wsum += n_edge * dc;
            hsum += n_edge * dc * E;
            used += n_edge;
        }
    }
    double d_bg = use_tilt ? s.d_of_c(1) : 1;
    double rho_sel = wsum / (wsum + std::max(m - used, 0.0) * d_bg);
    double H_rep = wsum > 0 ? hsum / wsum : 0;
    return (1 - rho_sel) * H_coll(n / std::pow(4.0, k) * coll_agree) + rho_sel * H_rep;
}

Prediction predict(const Scheme& s, const Families& fam, double m, double n, double theta, bool outer_survival, bool edges, bool shared_mutations) {
    double b = s.b;
    auto b_fn = s.b_of_L ? s.b_of_L : [b](double) { return b; };
    auto one = [](double) { return 1.0; };
    auto H = [&](double C, std::function<double(double)> bf, bool tilt, double ca) {
        return expected_H(s, fam, m, n, theta, outer_survival, edges, shared_mutations, C, bf, tilt, ca);
    };
    double H_full = H(s.C, b_fn, true, s.coll_agree);
    double H_genome = H(1, one, false, 1), H_tilt = H(1, one, true, 1), H_noctx = H(s.C, one, true, s.coll_agree);
    return {H_full, H_genome, H_tilt - H_genome, H_tilt - H_noctx, H_noctx - H_full};
}

// Theorem 1 under a Zipf law fitted to a real genome, theta = 0: repeated positions have copy
// number c in [2, c_max] with weight c^(1 - alpha), tilted by d(c); copy lengths L (bases) have
// position weights `lengths`, a position lies in the boundary zone with probability
// beta(L) = min(1, 2r / (L - k + 1)); a fraction rho of positions is repeated.
double predict_zipf(const Scheme& s, double alpha, long c_max, const std::map<long, double>& lengths, double rho, double n) {
    double total = 0, beta = 0;
    for (auto& l : lengths) total += l.second;
    for (auto& l : lengths) beta += l.second / total * std::min(1.0, 2.0 * s.r / std::max(l.first - s.k + 1, 1L));
    double num = 0, den = 0, w_all = 0;
    for (long c = 2; c <= c_max; c++) {
        double w = std::pow((double)c, 1 - alpha), dc = s.d_of_c ? s.d_of_c(c) : 1;
        num += w * dc * ((1 - beta) * Phi(c - 1, s.C) + beta * Phi(c - 1, s.b));
        den += w * dc, w_all += w;
    }
    double d_rep = den / w_all, d_bg = s.d_of_c ? s.d_of_c(1) : 1;  // selection rate in repeats / background
    double rho_sel = rho * d_rep / (rho * d_rep + (1 - rho) * d_bg);
    return (1 - rho_sel) * H_coll(n / std::pow(4.0, s.k) * s.coll_agree) + rho_sel * num / den;
}

// ---------------------------------------------------------------------------
// Zipf closed forms (Section 7)
// ---------------------------------------------------------------------------
double zeta(double s, double q) {  // Euler-Maclaurin summation
    if (s == 1) return Inf;
    if (s < 1 || q <= 0) return NaN;
    const int N = 20;
    static const double B2j[] = {1.0 / 6, -1.0 / 30, 1.0 / 42, -1.0 / 30, 5.0 / 66,
                                 -691.0 / 2730, 7.0 / 6, -3617.0 / 510, 43867.0 / 798, -174611.0 / 330};
    double sum = 0, a = q + N;
    for (int n = 0; n < N; n++) sum += std::pow(q + n, -s);
    sum += std::pow(a, 1 - s) / (s - 1) + 0.5 * std::pow(a, -s);
    double fact = 1, poch = s, apow = std::pow(a, -s - 1);
    for (int j = 1; j <= 10; j++) {
        fact *= (2.0 * j - 1) * (2.0 * j);
        double term = B2j[j - 1] / fact * poch * apow;
        sum += term;
        if (std::fabs(term) < 1e-17 * std::fabs(sum)) break;
        poch *= (s + 2 * j - 1) * (s + 2 * j);
        apow /= a * a;
    }
    return sum;
}

static double zeta_prime(double s) { return (zeta(s + 1e-5) - zeta(s - 1e-5)) / 2e-5; }

double H_rep_closed(double alpha) {
    return alpha <= 2 ? Inf : -zeta_prime(alpha - 1) / (zeta(alpha - 1) * LN2);
}

double H_rep_c2(double alpha) { return (-zeta_prime(alpha - 1) / LN2) / (zeta(alpha - 1) - 1); }

double H_rep_truncated(double alpha, long c_max, long c_min) {
    double num = 0, den = 0;
    for (long c = c_min; c <= c_max; c++) {
        double P = std::pow((double)c, 1 - alpha);
        num += P * std::log2((double)c);
        den += P;
    }
    return num / den;
}

ClosedForms closed_forms(double alpha) {
    double z1 = zeta(alpha - 1);
    return {H_rep_closed(alpha), alpha > 2 ? H_rep_c2(alpha) : Inf, alpha > 2 ? 1 / z1 : 0,
            alpha > 1 ? 1 / zeta(alpha) : 0, alpha > 3 ? zeta(alpha - 2) / z1 : Inf};
}

// scipy.optimize.brentq, step for step
static double brentq(const std::function<double(double)>& f, double xa, double xb) {
    const double xtol = 2e-12, rtol = 8.881784197001252e-16;
    double xpre = xa, xcur = xb, xblk = 0, fblk = 0, spre = 0, scur = 0;
    double fpre = f(xpre), fcur = f(xcur);
    if (fpre == 0) return xpre;
    if (fcur == 0) return xcur;
    for (int i = 0; i < 100; i++) {
        if (fpre != 0 && fcur != 0 && std::signbit(fpre) != std::signbit(fcur)) {
            xblk = xpre, fblk = fpre;
            spre = scur = xcur - xpre;
        }
        if (std::fabs(fblk) < std::fabs(fcur)) {
            xpre = xcur, xcur = xblk, xblk = xpre;
            fpre = fcur, fcur = fblk, fblk = fpre;
        }
        double delta = (xtol + rtol * std::fabs(xcur)) / 2, sbis = (xblk - xcur) / 2;
        if (fcur == 0 || std::fabs(sbis) < delta) return xcur;
        if (std::fabs(spre) > delta && std::fabs(fcur) < std::fabs(fpre)) {
            double stry;
            if (xpre == xblk) {
                stry = -fcur * (xcur - xpre) / (fcur - fpre);
            } else {
                double dpre = (fpre - fcur) / (xpre - xcur), dblk = (fblk - fcur) / (xblk - xcur);
                stry = -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre));
            }
            if (2 * std::fabs(stry) < std::min(std::fabs(spre), 3 * std::fabs(sbis) - delta)) spre = scur, scur = stry;
            else spre = scur = sbis;
        } else {
            spre = scur = sbis;
        }
        xpre = xcur, fpre = fcur;
        xcur += std::fabs(scur) > delta ? scur : (sbis > 0 ? delta : -delta);
        fcur = f(xcur);
    }
    return xcur;
}

// scipy.optimize.minimize_scalar(method="bounded"), step for step
static double fminbound(const std::function<double(double)>& f, double a, double b) {
    const double xatol = 1e-5, sqrt_eps = std::sqrt(2.2e-16), golden = 0.5 * (3 - std::sqrt(5.0));
    double fulc = a + golden * (b - a), nfc = fulc, xf = fulc, rat = 0, e = 0, x = xf;
    double fx = f(x), ffulc = fx, fnfc = fx, xm = 0.5 * (a + b);
    double tol1 = sqrt_eps * std::fabs(xf) + xatol / 3, tol2 = 2 * tol1;
    auto sign = [](double v) { return v > 0 ? 1.0 : v < 0 ? -1.0 : 0.0; };
    for (int num = 1; std::fabs(xf - xm) > tol2 - 0.5 * (b - a) && num < 500; num++) {
        bool golden_step = true;
        if (std::fabs(e) > tol1) {  // parabolic fit
            golden_step = false;
            double r = (xf - nfc) * (fx - ffulc), q = (xf - fulc) * (fx - fnfc);
            double p = (xf - fulc) * q - (xf - nfc) * r;
            q = 2 * (q - r);
            if (q > 0) p = -p;
            q = std::fabs(q);
            r = e;
            e = rat;
            if (std::fabs(p) < std::fabs(0.5 * q * r) && p > q * (a - xf) && p < q * (b - xf)) {
                rat = (p + 0.0) / q;
                x = xf + rat;
                if (x - a < tol2 || b - x < tol2) rat = tol1 * (sign(xm - xf) + (xm - xf == 0));
            } else {
                golden_step = true;
            }
        }
        if (golden_step) {
            e = xf >= xm ? a - xf : b - xf;
            rat = golden * e;
        }
        x = xf + (sign(rat) + (rat == 0)) * std::max(std::fabs(rat), tol1);
        double fu = f(x);
        if (fu <= fx) {
            if (x >= xf) a = xf;
            else b = xf;
            fulc = nfc, ffulc = fnfc, nfc = xf, fnfc = fx, xf = x, fx = fu;
        } else {
            if (x < xf) a = x;
            else b = x;
            if (fu <= fnfc || nfc == xf) fulc = nfc, ffulc = fnfc, nfc = x, fnfc = fu;
            else if (fu <= ffulc || fulc == xf || fulc == nfc) fulc = x, ffulc = fu;
        }
        xm = 0.5 * (a + b);
        tol1 = sqrt_eps * std::fabs(xf) + xatol / 3;
        tol2 = 2 * tol1;
    }
    return xf;
}

double estimate_alpha(double H_rep, bool c_ge_2) {
    if (!std::isfinite(H_rep)) return NaN;
    auto g = [&](double a) { return (c_ge_2 ? H_rep_c2(a) : H_rep_closed(a)) - H_rep; };
    if (g(2.0005) < 0 || g(30) > 0) return NaN;  // outside what any alpha in (2, 30) produces
    return brentq(g, 2.0005, 30);
}

double fit_zipf_mle(const std::vector<double>& hist, long c_min, bool position_weighted, long c_max) {
    double n = 0, sl = 0;
    long c_obs = 0;
    for (long c = c_min; c < (long)hist.size(); c++)
        if (hist[c] > 0) n += hist[c], sl += hist[c] * std::log((double)c), c_obs = c;
    if (n < 10) return NaN;
    if (c_max == 0) c_max = c_obs;
    std::function<double(double)> obj;  // negative log-likelihood
    double lo;
    if (c_max < 0) {
        obj = [&](double a) { return n * std::log(zeta(a, (double)c_min)) + a * sl; };
        lo = 1.01;
    } else {
        obj = [&](double a) {
            double s = 0;
            for (long c = c_min; c <= c_max; c++) s += std::exp(-a * std::log((double)c));
            return n * std::log(s) + a * sl;
        };
        lo = -2;
    }
    return fminbound(obj, lo, 10) + (position_weighted ? 1 : 0);
}

// ---------------------------------------------------------------------------
// calibration (ue/calibrate.py): pairs of mutated copies between random flanks
// ---------------------------------------------------------------------------
Sketch calibration_sketch(const std::string& key, const Settings& s) {
    return [key, s](const std::string& seq) {
        Settings t = s;
        Counts own;
        if (key == "weighted-min") {
            std::vector<uint64_t> h;
            for (const Seed& x : all_kmers(seq, s.k)) h.push_back(x.hash);
            own = make_counts(h);
            t.genome = &own;
        }
        return select_seeds(key, seq, t, 0);
    };
}

static std::string random_bases(std::mt19937_64& rng, size_t n) {
    std::string s(n, 'A');
    for (auto& c : s) c = "ACGT"[rng() >> 62];
    return s;
}

static std::string mutate(std::mt19937_64& rng, std::string s, double theta) {
    const std::string acgt = "ACGT";
    for (auto& c : s)
        if ((rng() >> 11) * 0x1p-53 < theta) c = acgt[(acgt.find(c) + 1 + rng() % 3) % 4];
    return s;
}

// selected flag and seed identity at every position
static void positions(const std::string& s, Sketch select, Sketch space, std::vector<uint8_t>& sel, std::vector<uint64_t>& id) {
    sel.assign(s.size(), 0);
    id.assign(s.size(), 0);
    for (const Seed& x : space(s)) id[x.start] = x.hash;
    for (const Seed& x : select(s)) sel[x.start] = 1;
}

Calibration calibrate(Sketch select, Sketch space, int r, int span, int k, double theta, int L, int reps, uint64_t seed) {
    std::mt19937_64 rng(seed);
    int flank = 2 * (r + 60) + k;
    long agree_int[2] = {0, 0}, agree_bd[2] = {0, 0}, surv[2] = {0, 0}, dens[2] = {0, 0};
    std::map<int, long> prof_hit, prof_tot;
    bool has_interior = L - k + 1 > 2 * r + 10;
    auto enough = [&]() {  // C, b and survival each well determined
        long m = surv[1] / k;
        if (has_interior) m = std::min(m, agree_int[1]);
        if (r > 0) m = std::min(m, agree_bd[1]);
        return m >= 1000;
    };
    std::vector<uint8_t> selA, selB, sel0;
    std::vector<uint64_t> idA, idB, id0;
    for (int done = 0; done < reps || (done < 3000 && !enough()); done++) {
        std::string anc = random_bases(rng, L);
        std::string A = mutate(rng, anc, theta), B = mutate(rng, anc, theta);
        auto build = [&](const std::string& core) { return random_bases(rng, flank) + core + random_bases(rng, flank); };
        positions(build(A), select, space, selA, idA);
        positions(build(B), select, space, selB, idB);
        positions(build(anc), select, space, sel0, id0);
        for (int i = flank; i <= flank + L - span; i++) {  // seed fully inside the copy
            int dist = std::min(i - flank, flank + L - span - i);
            surv[0] += idA[i] == id0[i];
            surv[1]++;
            dens[0] += selA[i];
            dens[1]++;
            if (!selA[i] || idA[i] != idB[i]) continue;
            if (dist >= r) {
                agree_int[0] += selB[i], agree_int[1]++;
            } else {
                agree_bd[0] += selB[i], agree_bd[1]++;
                prof_hit[dist] += selB[i], prof_tot[dist]++;
            }
        }
    }
    Calibration cal;
    cal.p = (double)surv[0] / surv[1];
    cal.C = agree_int[1] ? (double)agree_int[0] / agree_int[1] : 1;
    cal.b = agree_bd[1] ? (double)agree_bd[0] / agree_bd[1] : 1;
    cal.d = (double)dens[0] / dens[1];
    for (auto& t : prof_tot) cal.b_profile[t.first] = (double)prof_hit[t.first] / t.second;
    cal.n_interior = agree_int[1];
    cal.n_boundary = agree_bd[1];
    return cal;
}

std::function<double(double)> b_of_length(Sketch select, Sketch space, int r, int span, int k, double theta, std::vector<int> lengths, int reps, uint64_t seed) {
    std::sort(lengths.begin(), lengths.end());
    lengths.erase(std::unique(lengths.begin(), lengths.end()), lengths.end());
    std::vector<double> xs, ys;
    for (int L : lengths) {
        xs.push_back(L);
        ys.push_back(calibrate(select, space, r, span, k, theta, L, std::max(2, reps * 600 / L), seed + L).b);
    }
    return [xs, ys](double L) {  // np.interp
        if (L <= xs.front()) return ys.front();
        if (L >= xs.back()) return ys.back();
        size_t j = std::upper_bound(xs.begin(), xs.end(), L) - xs.begin();
        return ys[j - 1] + (L - xs[j - 1]) / (xs[j] - xs[j - 1]) * (ys[j] - ys[j - 1]);
    };
}

std::vector<double> context_kernel(Sketch select, int r, int k, long hits, uint64_t seed) {
    std::mt19937_64 rng(seed);
    const std::string acgt = "ACGT";
    auto other = [&](char c) { return acgt[(acgt.find(c) + 1 + rng() % 3) % 4]; };
    auto selected = [&](const std::string& s, uint32_t pos) {
        for (const Seed& x : select(s))
            if (x.start == pos) return true;
        return false;
    };
    const int pad = r + 30;  // beyond every window that can select the seed
    std::vector<double> q;
    for (int e = 0; e < r; e++) {
        long sel = 0, agree = 0;
        for (long it = 0; sel < hits && it < 1000 * hits; it++) {
            std::string core = random_bases(rng, k + 2 * e);
            std::string A = random_bases(rng, pad) + core + random_bases(rng, pad);
            std::string B = random_bases(rng, pad) + core + random_bases(rng, pad);
            B[pad - 1] = other(A[pad - 1]), B[pad + k + 2 * e] = other(A[pad + k + 2 * e]);  // differ next to it
            if (!selected(A, pad + e)) continue;
            sel++, agree += selected(B, pad + e);
        }
        q.push_back(sel ? (double)agree / sel : NaN);
    }
    return q;
}

// random unit of minimal period exactly p
static std::string periodic_unit(std::mt19937_64& rng, int p) {
    for (;;) {
        std::string u = random_bases(rng, p);
        bool shorter = false;
        for (int q = 1; q < p && !shorter; q++)
            if (p % q == 0) {
                bool per = true;
                for (int j = q; j < p && per; j++) per = u[j] == u[j - q];
                shorter = per;
            }
        if (!shorter) return u;
    }
}

TandemRates tandem_rates(SketchHS select, int r, int k, int p_max, int reps_per_length, uint64_t seed) {
    std::mt19937_64 rng(seed);
    TandemRates t;
    t.r = r, t.p_max = p_max = std::min(p_max, k / 2);
    size_t cells = (size_t)(p_max + 1) * (r + 1) * (r + 1);
    std::vector<double> hit(cells, 0), tot(cells, 0);
    double bg_hit = 0, bg_tot = 0;
    const int flank = k + r + 30;
    std::vector<int> lengths;  // k-mer positions per run: every two-edge case, then a long run
    for (int l = 1; l <= 2 * r + 2; l++) lengths.push_back(l);
    lengths.push_back(2 * r + 40);
    uint64_t hs = 0;
    std::vector<uint8_t> sel;
    for (int p = 1; p <= p_max; p++)
        for (int l : lengths)
            for (int rep = 0; rep < reps_per_length; rep++, hs++) {
                std::string u = periodic_unit(rng, p), run;
                int phase = (int)(rng() % p);
                for (int j = 0; j < l + k - 1; j++) run += u[(j + phase) % p];
                std::string s = random_bases(rng, flank) + run + random_bases(rng, flank);
                sel.assign(s.size(), 0);
                for (const Seed& x : select(s, hs)) sel[x.start] = 1;
                // the stretch of periodic k-mers that contains the run (flanks may extend it by chance)
                std::vector<PeriodicKmer> pk = periodic_kmers(s, k, p_max);
                size_t a = 0;
                while (a < pk.size() && pk[a].start < (uint32_t)flank) a++;
                if (a == pk.size() || pk[a].start != (uint32_t)flank) continue;
                size_t lo = a, hi = a;
                while (lo > 0 && pk[lo - 1].start + 1 == pk[lo].start && pk[lo - 1].p == pk[a].p) lo--;
                while (hi + 1 < pk.size() && pk[hi + 1].start == pk[hi].start + 1 && pk[hi + 1].p == pk[a].p) hi++;
                for (size_t j = lo; j <= hi; j++) {
                    int dl = std::min<int>(j - lo, r), dr = std::min<int>(hi - j, r);
                    size_t i = ((size_t)pk[j].p * (r + 1) + dl) * (r + 1) + dr;
                    hit[i] += sel[pk[j].start], tot[i]++;
                }
                // aperiodic flanks beyond every window that reaches the run
                for (int i = r; i + k + r < flank; i++) bg_hit += sel[i], bg_tot++;
                for (size_t i = flank + run.size() + r; i + k + r <= s.size(); i++) bg_hit += sel[i], bg_tot++;
            }
    t.base = bg_tot ? bg_hit / bg_tot : 0;
    t.rate.assign(cells, NaN), t.n = tot;
    for (size_t i = 0; i < cells; i++)
        if (tot[i] > 0) t.rate[i] = hit[i] / tot[i];
    for (int p = 1; p <= p_max; p++)  // unsampled cells (rare): the mirror cell, else the interior
        for (int dl = 0; dl <= r; dl++)
            for (int dr = 0; dr <= r; dr++) {
                size_t i = ((size_t)p * (r + 1) + dl) * (r + 1) + dr, m = ((size_t)p * (r + 1) + dr) * (r + 1) + dl,
                       c = ((size_t)p * (r + 1) + r) * (r + 1) + r;
                if (std::isnan(t.rate[i])) t.rate[i] = !std::isnan(t.rate[m]) ? t.rate[m] : t.rate[c];
            }
    return t;
}

double collision_agreement(Sketch select, int r, int k, int reps, uint64_t seed) {
    if (r == 0) return 1;
    std::mt19937_64 rng(seed);
    int pad = r + 5;
    long hit = 0, tot = 0;
    std::vector<uint8_t> sel;
    std::vector<uint64_t> id;
    auto selected = [&](const std::string& core) {
        std::string s = random_bases(rng, pad) + core + random_bases(rng, pad);
        for (const Seed& x : select(s))
            if (x.start == (uint32_t)pad) return true;
        return false;
    };
    for (int i = 0; i < reps; i++) {
        std::string core = random_bases(rng, k);
        if (selected(core)) tot++, hit += selected(core);
    }
    return tot ? (double)hit / tot : NaN;
}
