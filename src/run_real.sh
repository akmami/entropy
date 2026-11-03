#!/usr/bin/env bash
# Real genomes: measure every scheme, predict H from the genome's repeat structure, compare.
#
#   src/run_real.sh
#
# Every run starts from scratch: the outputs of the previous run are deleted (the genomes are kept).
#   1. download every genome missing from genomes/real/<name>.fa.gz
#   2. for each genome:
#        src/copies   suffix array + LCP -> genomes/real/<name>.copies.tsv
#        src/ue real  every scheme, predictions, alpha, k-sweep -> results/real/<name>.json (log: <name>.log)
#   3. src/plot.py -> results/real/ (tables and figure)
#
# Human (GRCh38): copies needs about 80 GB of memory and 1-2 hours, ue about 150 GB and 6 hours.
set -eu
cd "$(dirname "$0")/.."

K=21                    # k-mer length
W=11                    # minimizer window
KS=15,17,19,21,25,31    # k values of the k-sweep
DIR=genomes/real
OUT=results/real

GENOMES="
ecoli  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
yeast  https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
dmel   https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
poplar https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/775/GCF_000002775.5_P.trichocarpa_v4.1/GCF_000002775.5_P.trichocarpa_v4.1_genomic.fna.gz
human  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz
"

rm -rf "$OUT" "$DIR"/*.copies.tsv
mkdir -p "$DIR" "$OUT"

make -C src
RADII=$(src/ue settings --k $K --w $W | awk '$1 == "radii" {print $2}')  # context radii copies must count

# 1. download: a genome goes into the list if it is already there or downloads successfully
genomes=()
while read -r name url; do
    [ -n "$name" ] || continue
    fa=$DIR/$name.fa.gz
    if [ ! -f "$fa" ]; then
        echo ">> $name: downloading $url"
        if ! curl -fL --retry 5 -o "$fa.part" "$url"; then
            rm -f "$fa.part"
            echo "!! $name: download failed, skipped" >&2
            continue
        fi
        mv "$fa.part" "$fa"
    fi
    genomes+=("$name")
done <<< "$GENOMES"

# 2. process
process() {
    local name=$1 fa=$DIR/$1.fa.gz log=$OUT/$1.log
    echo ">> $name"
    src/copies "$fa" $K $RADII $KS > "$log" 2>&1 &&
        src/ue real "$fa" --k $K --w $W --ks $KS --out "$OUT/$name.json" >> "$log" 2>&1 ||
        echo "!! $name failed, see $log" >&2
}

for name in "${genomes[@]}"; do
    process "$name"
done

# 3. plot
python3 src/plot.py
echo "done: $OUT/"
