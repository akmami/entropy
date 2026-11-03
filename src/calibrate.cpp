// calibrate: the calibrated quantities of Theorem 1 for every context-dependent scheme, from pairs of
// simulated repeat copies (no genome). One file for all experiments: ue panel reads it (--calibration).
//
//   src/calibrate --ks "21 31" --ws "11 41" --thetas "0 0.05" --out results/simulated/calibrate.tsv
//
// One row per scheme, k, w and theta (lists take spaces or commas). The schemes are the ones with a
// context, the minimizers by default; context-free schemes have C = b = 1. A copy pair is two copies
// of a random template, each with substitutions at rate theta, between random flanks
// (theory.cpp: calibrate()). Columns:
//   r              context radius (w - 1; the mod-minimizer: its own window - 1)
//   d              density on random sequence
//   p              seed survival: the seed at a position of a copy equals the template's
//   C              co-selection inside a copy (>= r from both edges):
//                  P[selected in copy B | selected in copy A, same seed]
//   b              the same within r of an edge (boundary agreement), copies of COPY_LENGTH
//   collision      P[a chance copy of a selected k-mer, between random flanks, is selected too]
//   n_interior, n_boundary   selected pairs behind C and b
//   b_by_length    L:b for short copies, whose positions can be near both edges; longer copies have b
//   b_by_distance  distance:b, b by distance 0, 1, ..., r-1 to the nearest edge (the boundary kernel)
//
// build: make -C src calibrate
#include <chrono>
#include <cstdio>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "sketches.h"
#include "theory.h"

static const int COPY_LENGTH = 2000;  // copies for C and b
static const int REPS = 20;           // least number of copy pairs; more until C and b are well determined

static const char* USAGE =
    "usage: calibrate --ks K... --ws W... --thetas THETA... --out FILE [--schemes A,B,...]\n"
    "       [schemes: random-min,lex-min,weighted-min,mod-min]\n";

// "21 11,13" -> {"21", "11", "13"}
static std::vector<std::string> list(std::string s) {
    for (char& c : s)
        if (c == ',') c = ' ';
    std::istringstream in(s);
    std::vector<std::string> out;
    for (std::string x; in >> x;) out.push_back(x);
    return out;
}

static std::string fmt(double x) {
    char b[32];
    snprintf(b, sizeof b, "%.6g", x);
    return b;
}

static double seconds() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

int main(int argc, char** argv) {
    std::map<std::string, std::string> opt = {{"--schemes", "random-min,lex-min,weighted-min,mod-min"}};
    const std::set<std::string> known = {"--ks", "--ws", "--thetas", "--out", "--schemes"};
    for (int i = 1; i < argc; i += 2) {
        if (!known.count(argv[i]) || i + 1 == argc) return fputs(USAGE, stderr), 2;
        opt[argv[i]] = argv[i + 1];
    }
    for (const char* o : {"--ks", "--ws", "--thetas", "--out"})
        if (!opt.count(o)) return fputs(USAGE, stderr), 2;

    std::string out = opt["--out"];
    FILE* f = fopen((out + ".tmp").c_str(), "w");
    if (!f) return fprintf(stderr, "calibrate: cannot write %s\n", out.c_str()), 2;
    fprintf(f, "scheme\tk\tw\ttheta\tr\td\tp\tC\tb\tcollision\tn_interior\tn_boundary\tb_by_length\tb_by_distance\n");
    for (auto& k_str : list(opt["--ks"]))
        for (auto& w_str : list(opt["--ws"])) {
            Settings s;
            s.k = std::stoi(k_str), s.w = std::stoi(w_str), s.d = 2.0 / (s.w + 1);
            s.modmin_w = modmin_window(s.k, s.d);
            for (auto& theta_str : list(opt["--thetas"]))
                for (auto& key : list(opt["--schemes"])) {
                    int r = radius(key, s), sn = span(key, s);
                    if (r == 0) continue;  // context-free: C = b = 1
                    double theta = std::stod(theta_str), t0 = seconds();
                    Sketch sel = calibration_sketch(key, s);
                    Sketch sp = [key, s](const std::string& seq) { return seed_space(key, seq, s, 0); };

                    Calibration c = calibrate(sel, sp, r, sn, s.k, theta, COPY_LENGTH, REPS);
                    double coll = collision_agreement(sel, r, s.k, 3000);
                    std::vector<int> lengths;  // from the shortest copy with a whole window to 6 (r + 1) + 2k
                    int lo = r + s.k, hi = 6 * (r + 1) + 2 * s.k;
                    for (int i = 0; i < 6; i++) lengths.push_back(lo + (hi - lo) * i / 5);
                    auto b_of_L = b_of_length(sel, sp, r, sn, s.k, theta, lengths, 8);

                    std::string by_length, by_distance;
                    for (int L : lengths) by_length += (by_length.empty() ? "" : ",") + std::to_string(L) + ":" + fmt(b_of_L(L));
                    for (auto& kv : c.b_profile)
                        by_distance += (by_distance.empty() ? "" : ",") + std::to_string(kv.first) + ":" + fmt(kv.second);
                    fprintf(f, "%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%ld\t%ld\t%s\t%s\n", key.c_str(), s.k, s.w,
                            fmt(theta).c_str(), r, fmt(c.d).c_str(), fmt(c.p).c_str(), fmt(c.C).c_str(), fmt(c.b).c_str(),
                            fmt(coll).c_str(), c.n_interior, c.n_boundary, by_length.c_str(), by_distance.c_str());
                    fprintf(stderr, "%-12s k=%-2d w=%-2d theta=%-4s C=%.3f b=%.3f (%.1fs)\n", key.c_str(), s.k, s.w,
                            fmt(theta).c_str(), c.C, c.b, seconds() - t0);
                }
        }
    fclose(f);
    rename((out + ".tmp").c_str(), out.c_str());
    fprintf(stderr, "wrote %s\n", out.c_str());
    return 0;
}
