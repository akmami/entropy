# Uniqueness entropy of k-mer sketches

Code for the paper *A Scheme-Comparable Formulation of Uniqueness Entropy for k-mer Sketches on Repetitive Genomes*. It measures the uniqueness entropy of k-mer selection schemes (minimizers, mod-minimizers, syncmers, FracMinHash and others) on simulated and real genomes, and compares it with the prediction of Theorem 1.

## Layout

| path | contents |
|---|---|
| `src/` | C++ tools (`ue`, `calibrate`, `copies`), the genome simulator, plotting and run scripts |
| `manuscript/` | the paper |
| `results/` | outputs, written by the run scripts |
| `misc/` | an earlier Python implementation, kept for reference |

## Build

Needs a C++17 compiler and zlib; `copies` also needs cmake (it downloads and builds libdivsufsort on the first `make`). The Python scripts need the packages in `requirements.txt`.

```bash
make -C src
```

## Run

```bash
src/run_sim.sh    # simulated genomes  -> results/simulated/
src/run_real.sh   # real genomes       -> results/real/
```

Each script runs everything from scratch and ends with `src/plot.py`, which writes the tables and figures. Settings are at the top of each script. The simulations take about 3.5 CPU hours; the human genome needs about 150 GB of memory.
