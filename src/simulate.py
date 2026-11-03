#!/usr/bin/env python3
"""Planted-Zipf genomes (Definition: planted Zipf repeats), --reps genomes of one setting.

    python3 src/simulate.py --out genomes/simulated --n 20000000 --alpha 2.5 --rho 0.05 \
        --theta 0 --Lmin 101 --Lmax 500 --cmax 2000 --reps 10

Every argument is required; src/run_sim.sh holds the settings. Each genome is named by its
setting, n<Mb>M_a<alpha>_t<theta>_r<replicate>, and written as two files:

    n20M_a2.5_t0_r0.fa          the genome, one record
    n20M_a2.5_t0_r0.truth.tsv   '#key<TAB>value' lines with the parameters, then one row per
                                planted copy, by position: family, copies, length, start (0-based)

The name also sets the random seed (its CRC32), so a name always gives the same genome. A genome
that already exists is kept, so a setting shared by several sweeps is made once. The path of
every genome is printed.

The model: a uniform random background into which repeat families are planted until their
copies cover rho * n bases. A family has c copies, P[c] ~ c^-alpha on [2, cmax], of length L
uniform in [Lmin, Lmax]. Each copy is the family's random template with every base substituted
with probability theta. Copies go in random order, with at least GAP background bases between them.
"""
import argparse
import os
import zlib

import numpy as np

GAP = 64  # least background between two copies


def planted_genome(rng, n, alpha, rho, theta, Lmin, Lmax, cmax):
    """The genome (uint8, 0..3 = ACGT) and its copies as rows (family, copies, length, start)."""
    # families, until their copies fill the repeat budget; the last one is cut to fit
    c_values = np.arange(2, cmax + 1)
    p = c_values.astype(float) ** -alpha
    p /= p.sum()
    budget, used, families = int(rho * n), 0, []
    while True:
        L = int(rng.integers(Lmin, Lmax + 1))
        c = min(int(rng.choice(c_values, p=p)), (budget - used) // L)
        if c < 2:
            break
        families.append((c, rng.integers(0, 4, size=L, dtype=np.uint8)))
        used += c * L

    # copies: the template with substitutions, each to one of the other three bases
    copies = []
    for j, (c, template) in enumerate(families):
        for _ in range(c):
            copy = template.copy()
            hit = rng.random(len(copy)) < theta
            copy[hit] = (copy[hit] + rng.integers(1, 4, size=hit.sum())) % 4
            copies.append((j, copy))

    # background gaps: GAP each, plus the remaining background split at random
    free = n - used - GAP * (len(copies) + 1)
    assert free >= 0, "the copies do not fit into the genome"
    cuts = np.sort(rng.integers(0, free + 1, size=len(copies)))
    gaps = np.diff(np.concatenate(([0], cuts, [free]))) + GAP

    seq = rng.integers(0, 4, size=n, dtype=np.uint8)  # background
    rows, pos = [], 0
    for gap, i in zip(gaps, rng.permutation(len(copies))):
        j, copy = copies[i]
        pos += gap
        seq[pos : pos + len(copy)] = copy
        rows.append((j, families[j][0], len(copy), pos))
        pos += len(copy)
    return seq, rows


def write_fasta(path, seq):
    text = np.frombuffer(b"ACGT", dtype=np.uint8)[seq].tobytes()
    with open(path, "wb") as f:
        f.write(b">planted\n")
        for i in range(0, len(text), 80):
            f.write(text[i : i + 80] + b"\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="output folder")
    ap.add_argument("--n", type=int, required=True, help="genome length")
    ap.add_argument("--alpha", type=float, required=True, help="Zipf exponent of the copy numbers")
    ap.add_argument("--rho", type=float, required=True, help="fraction of bases in repeat copies")
    ap.add_argument("--theta", type=float, required=True, help="substitution rate of a copy")
    ap.add_argument("--Lmin", type=int, required=True, help="shortest copy")
    ap.add_argument("--Lmax", type=int, required=True, help="longest copy")
    ap.add_argument("--cmax", type=int, required=True, help="largest copy number")
    ap.add_argument("--reps", type=int, required=True, help="genomes of this setting")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    for rep in range(a.reps):
        name = f"n{a.n / 1e6:g}M_a{a.alpha:g}_t{a.theta:g}_r{rep}"
        path = os.path.join(a.out, name)
        if not os.path.exists(path + ".fa"):
            seed = zlib.crc32(name.encode())
            seq, rows = planted_genome(np.random.default_rng(seed), a.n, a.alpha, a.rho, a.theta,
                                       a.Lmin, a.Lmax, a.cmax)
            params = dict(n=a.n, alpha=a.alpha, rho=a.rho, theta=a.theta, L_min=a.Lmin, L_max=a.Lmax,
                          c_max=a.cmax, rep=rep, seed=seed, families=len({r[0] for r in rows}),
                          rho_realized=sum(r[2] for r in rows) / a.n)
            with open(path + ".truth.tsv", "w") as f:
                for key, value in params.items():
                    f.write(f"#{key}\t{value}\n")
                f.write("family\tcopies\tlength\tstart\n")
                for row in rows:
                    f.write("\t".join(map(str, row)) + "\n")
            write_fasta(path + ".fa", seq)  # last: an existing .fa means the genome is complete
        print(path + ".fa")


if __name__ == "__main__":
    main()
