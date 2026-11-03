#!/usr/bin/env bash
# Simulated genomes: measure every scheme and predict H from the planted repeats.
#
#   src/run_sim.sh
#
# Every run starts from scratch: genomes/simulated/ and results/simulated/ are deleted.
#   1. src/simulate.py  10 genomes per setting -> genomes/simulated/<name>.fa and <name>.truth.tsv
#                       base: n = 20 Mb, alpha = 2.5, rho = 0.05, theta = 0, copies of 101-500 bp, c_max = 2000
#                       sweeps, each changing one base value: alpha, theta, n
#   2. src/calibrate    C and b of every context-dependent scheme at every k, w and theta, from simulated
#                       copy pairs -> results/simulated/calibrate.tsv, which every ue panel run reads
#   3. src/ue panel     every genome at k = 21, w = 11; the base genomes also at the other w and k
#                       -> results/simulated/<name>_k<k>_w<w>.json
#   4. src/ue           Zipf closed forms and alpha from a sketch, on the alpha sweep
#   5. src/plot.py      -> results/simulated/: summary.csv, predicted_vs_measured (one row per sweep:
#                       measured against predicted H and Delta), schemes (two schemes per row), schemes.csv
#
# THREADS runs at once (xargs -P), each on one core with about 2 GB per 20 Mb of genome. In all about
# 3.5 CPU hours: about 20 minutes with 32 threads.
set -eu
cd "$(dirname "$0")/.."

THREADS=10                                        # runs at once
K=21                                              # k-mer length
W=11                                              # minimizer window
WS="5 21 41 81"                                   # other w, with k = 21 (copies are >= w + k - 1 long)
KS="11 13 15 17 25 31"                            # other k, with w = 11
ALPHAS="1.5 2 2.5 3 3.5"                          # Zipf exponent of the copy numbers
THETAS="0 0.01 0.02 0.05 0.1"                     # substitution rate of a copy
NS="2500000 5000000 10000000 20000000 40000000"   # genome length
DIR=genomes/simulated
OUT=results/simulated

rm -rf "$DIR" "$OUT"
mkdir -p "$DIR" "$OUT"
make -C src ue calibrate

# 1. simulate: the base setting, changed by the arguments (argparse keeps an option's last value);
#    prints the genomes' paths
simulate() {
    python3 src/simulate.py --out $DIR --reps 10 --n 20000000 --alpha 2.5 --rho 0.05 --theta 0 \
        --Lmin 101 --Lmax 500 --cmax 2000 "$@"
}
base=$(simulate)
alpha_sweep=$(for a in $ALPHAS; do simulate --alpha $a; done)
for t in $THETAS; do simulate --theta $t; done
for n in $NS; do simulate --n $n; done

# 2. calibrate: one process per k, then one file (the header once)
for k in $K $KS; do echo $k; done |
    xargs -P $THREADS -I{} src/calibrate --ks {} --ws "$W $WS" --thetas "$THETAS" --out $OUT/calibrate_k{}.tsv
awk 'NR == 1 || FNR > 1' $OUT/calibrate_k*.tsv > $OUT/calibrate.tsv
rm $OUT/calibrate_k*.tsv

# 3. panel genome.fa k w, one line per run, THREADS at once
panel() {
    local name=$(basename $1 .fa)
    echo ">> $name k=$2 w=$3"
    src/ue panel $1 --k $2 --w $3 --hash-seeds 0,1,2 --calibration $OUT/calibrate.tsv --quiet --out $OUT/${name}_k$2_w$3.json
}
export -f panel
export OUT
{
    for fa in $DIR/*.fa; do 
        echo $fa $K $W; 
    done
    for fa in $base; do
        for w in $WS; do echo $fa $K $w; done
        for k in $KS; do echo $fa $k $W; done
    done
} | xargs -P $THREADS -L 1 bash -c 'panel "$@"' _

# 4. closed forms and alpha
src/ue closed-forms $alpha_sweep --k $K --quiet --out $OUT/closed_forms.json
src/ue estimate-alpha $alpha_sweep --k $K --quiet --out $OUT/estimate_alpha.json

# 5. plot
python3 src/plot.py
echo "done: $OUT/"
