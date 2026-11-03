#!/usr/bin/env python3
"""Tables and figures from the outputs of src/run_sim.sh and src/run_real.sh.

    python3 src/plot.py

No options: it redoes every part whose outputs exist.
  simulated genomes  results/simulated/<genome>_k<k>_w<w>.json -> results/simulated/
                       raw.csv                   one row per genome x k x w x hash seed x scheme
                       summary.csv               mean and sd over genomes and hash seeds, per setting x scheme
                       predicted_vs_measured     Theorem 1 against simulation, one row per sweep (alpha, theta,
                                                 n, w, k): measured against predicted H and Delta
                       schemes                   the same, one row per scheme, every setting of every sweep
                       schemes.csv               per scheme: values at the base setting, prediction errors
                       schemes.tex               the same as a LaTeX table (booktabs), for the paper
                       closed_forms.csv, estimate_alpha.csv
  real genomes       results/real/*.json -> results/real/  raw.csv, summary.csv, ue's tables, figure
"""
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter, LogLocator, MaxNLocator, NullFormatter  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# colors follow the scheme, never its rank; all k-mers and the oracle are references
COLORS = {"random-min": "#2a78d6", "lex-min": "#eb6834", "mod-min": "#1baf7a", "open-sync": "#eda100",
          "closed-sync": "#e87ba4", "fracminhash": "#008300", "weighted-min": "#4a3aa7", "randstrobe": "#e34948",
          "modmer": "#008300", "all": "#52514e", "oracle": "#0b0b0b"}
LABELS = {"all": "all k-mers", "random-min": "random minimizer", "lex-min": "lex minimizer",
          "weighted-min": "weighted minimizer", "mod-min": "mod-minimizer", "open-sync": "open syncmer",
          "closed-sync": "closed syncmer", "fracminhash": "FracMinHash", "modmer": "modmer",
          "randstrobe": "randstrobe", "oracle": "oracle"}
# the columns summary.csv averages (those a table has)
METRICS = ["d", "H", "H_true", "Delta", "U_N", "U_D", "Ehits", "H_min", "spurious_mass",
           "singleton_precision", "H_pred", "H_pred_crude", "H_genome_pred", "Delta_ctx_pred", "Delta_cons_pred",
           "Delta_order_pred", "H_bg", "H_rep", "frac_rep", "C_cal", "b_cal", "H_pred_thm", "H_pred_zipf", "H_pred_measured_d", "H_true_pred", "Delta_pred",
           "H_pred_no_tandem", "H_true_pred_no_tandem", "H_pred_zipf_distinct", "alpha_fit_distinct",
           "beta_emp", "beta_zipf", "agree_pred", "alpha_fit", "rho_fit",
           "unique_anchor_rate", "false_unique_rate", "no_unique_rate", "anchors_per_read", "hits_per_read",
           "hits_per_anchor"]
# the three predictions for real genomes
# real genomes, smallest first, and their names in figures
GENOMES = {"ecoli": "E. coli", "yeast": "yeast", "dmel": "fly", "poplar": "poplar", "human": "human"}
PRED = {"H_pred": "plug-in (copies and contexts)", "H_pred_thm": "Theorem 1, genome spectrum",
        "H_pred_zipf": "Theorem 1, fitted Zipf law"}

GRID = "#e4e3df"
plt.rcParams.update({"axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 11,
                     "axes.titlelocation": "left", "legend.frameon": False, "legend.fontsize": 8})


def load(pattern, table="rows"):
    """One table of every JSON file matching the pattern."""
    rows = []
    for path in sorted(glob.glob(pattern)):
        with open(path) as fh:
            rows += json.load(fh).get(table, [])
    return pd.DataFrame(rows)


def summarise(df, keys):
    """Mean and sd of every metric per (keys, scheme)."""
    g = df.groupby(keys + ["scheme"], sort=False)[[c for c in METRICS if c in df]]
    summ = pd.concat([g.mean().add_suffix("_mean"), g.std().add_suffix("_sd")], axis=1).reset_index()
    summ["n_runs"] = g.size().values
    return summ


def outdir(name):
    print(f">> results/{name}/")
    os.makedirs(f"results/{name}", exist_ok=True)
    return f"results/{name}"


def schemes_in(df):
    return [s for s in COLORS if s in set(df["scheme"])]


def save(fig, prefix):
    for ext in ("png", "pdf"):
        fig.savefig(f"{prefix}.{ext}", dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# simulated genomes
# ---------------------------------------------------------------------------
SIM = "results/simulated"
# the swept parameters, each with its axis label; a setting is one value of each
SWEEPS = [("alpha", "Zipf exponent α"), ("theta", "divergence θ"), ("n", "genome length n (Mb)"), ("w", "window w"),
          ("k", "k")]
SETTING = [p for p, _ in SWEEPS]
# the measured quantities and their Theorem 1 predictions: column, axis label, title, prediction
QUANTITIES = [("H", "H (bits)", "uniqueness entropy H", "H_pred"), ("Delta", "Δ (bits)", "false uniqueness Δ", "Delta_pred")]
# larger text for the figures, which are shrunk to the page width
# the real-genome figure is one row of 4 panels, shrunk more at page width: larger text still
WIDE_FONTS = {"font.size": 15, "axes.labelsize": 15, "axes.titlesize": 15, "xtick.labelsize": 14, "ytick.labelsize": 14,
              "legend.fontsize": 14, "legend.title_fontsize": 15}
PANEL_LABEL = dict(fontsize=17, fontweight="bold")  # (a), (b), ...: about 8 pt once the figure is at page width
FIG_FONTS = {"font.size": 13, "axes.labelsize": 13, "axes.titlesize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
             "legend.fontsize": 12, "legend.title_fontsize": 13, "figure.titlesize": 15}
# per-scheme figure: color = the sweep (categorical slots 1-5 in order), marker = its value, in order of the value;
# the base setting, which every sweep shares, is a gray star
SWEEP_COLOR = {"alpha": "#2a78d6", "theta": "#eb6834", "n": "#1baf7a", "w": "#eda100", "k": "#e87ba4"}
SWEEP_NAME = {"alpha": "α", "theta": "θ", "n": "n", "w": "w", "k": "k"}
BASE_COLOR = "#52514e"


def simulated():
    df = load(f"{SIM}/*_k*_w*.json")
    if df.empty:
        return
    out = outdir("simulated")
    df["Delta_pred"] = df["Delta_ctx_pred"] + df["Delta_cons_pred"]  # Corollary 1: boundary + interior loss
    df.to_csv(f"{out}/raw.csv", index=False)
    summ = summarise(df, SETTING)
    # the base setting: the value of each parameter that all the other sweeps keep, i.e. its most common one
    base = {p: summ[p].mode()[0] for p in SETTING}
    differs = [[p for p in SETTING if not np.isclose(r[p], base[p])] for _, r in summ.iterrows()]
    summ["sweep"] = [d[0] if d else "base" for d in differs]  # every setting changes at most one parameter
    summ.to_csv(f"{out}/summary.csv", index=False)
    with plt.rc_context(FIG_FONTS):
        predicted_vs_measured(summ, out)
        schemes(summ, out)
    scheme_table(summ, out)
    scheme_latex(summ, out)
    for table in ("closed_forms", "estimate_alpha"):
        load(f"{SIM}/{table}.json").to_csv(f"{out}/{table}.csv", index=False)


def predicted_vs_measured(summ, out):
    """Theorem 1 against simulation, a pair of panels per sweep (two sweeps per row, the last one centered below):
    measured (x) and predicted (y) H and Delta, one point per setting and scheme (mean over genomes and hash seeds);
    color: scheme, marker: the swept value, listed in one line under the pair."""
    d = summ[summ["H_pred_mean"].notna()]
    fig = plt.figure(figsize=(13, 11.5), layout="constrained")
    rows = fig.subfigures(3, 1, hspace=0.03)
    slots = [*rows[0].subfigures(1, 2, wspace=0.15), *rows[1].subfigures(1, 2, wspace=0.15),
             rows[2].subfigures(1, 3, width_ratios=[1, 2, 1])[1]]
    for letter, sf, (p, name) in zip("abcde", slots, SWEEPS):
        sf.suptitle(f"({letter})", x=0.01, ha="left", **PANEL_LABEL)
        sub = d[d["sweep"].isin([p, "base"])]
        marker = dict(zip(sorted(sub[p].unique()), "osD^vPX"))
        for ax, (col, label, title, pred) in zip(sf.subplots(1, len(QUANTITIES)), QUANTITIES):
            for _, r in sub.iterrows():
                ax.plot(r[f"{col}_mean"], r[f"{pred}_mean"], marker[r[p]], color=COLORS[r["scheme"]], ms=7,
                        mec="white", mew=0.5)
            v = sub[[f"{col}_mean", f"{pred}_mean"]]
            lo, hi = v.min().min(), v.max().max()
            pad = 0.05 * (hi - lo) if hi > lo else 0.01
            lim = (lo - pad, hi + pad)
            ax.plot(lim, lim, color="#8a8985", lw=1, ls="--", zorder=0)  # y = x
            ax.set(xlim=lim, ylim=lim, xlabel=f"measured {label}", ylabel=f"predicted {label}")
            linear_ticks(ax)
        # the sweep's name, then its values with their markers, in one line (n in Mb, as the name says)
        values = [plt.Line2D([], [], ls="none", label=name)]
        values += [plt.Line2D([], [], color="#0b0b0b", marker=m, mfc="none", ls="none", ms=8,
                              label=f"{v / 1e6:g}" if p == "n" else f"{v:g}") for v, m in marker.items()]
        legend = sf.legend(handles=values, loc="outside lower center", ncol=len(values), columnspacing=1.0,
                           handletextpad=0.3)
        legend.get_texts()[0].set_fontweight("bold")
    schemes = [plt.Line2D([], [], color=COLORS[s], marker="o", ls="none", ms=8, label=LABELS[s]) for s in schemes_in(d)]
    fig.legend(handles=schemes, loc="outside lower center", ncol=4)
    save(fig, f"{out}/predicted_vs_measured")


def linear_ticks(ax):
    """At most 4 ticks on x and 5 on y, written as plain numbers (0.0005, not 0.0005000)."""
    plain = FuncFormatter(lambda v, _: f"{v:g}")
    ax.xaxis.set_major_locator(MaxNLocator(3))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_formatter(plain)
    ax.yaxis.set_major_formatter(plain)


def log_axes(ax, decades):
    """Log x and y, ticks at powers of ten (and 2 and 5 times them over at most two decades; every other power over
    more than four), as plain numbers."""
    ax.set(xscale="log", yscale="log")
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(LogLocator(subs=(1, 2, 5)) if decades <= 2 else LogLocator(subs=(1,), numticks=3))
        axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        axis.set_minor_formatter(NullFormatter())


def schemes(summ, out):
    """Theorem 1 against simulation per scheme, two schemes side by side: measured (x) and predicted (y) H and Delta
    at every setting of every sweep (mean over genomes and hash seeds); color: the sweep, marker: its value, star: the
    base setting. H, and Delta where it is positive, on log axes."""
    d = summ[summ["H_pred_mean"].notna()]
    names = schemes_in(d)
    base = d[d["sweep"] == "base"].iloc[0]
    values = {p: sorted(d.loc[d["sweep"].isin([p, "base"]), p].unique()) for p in SWEEP_COLOR}
    markers = {p: dict(zip([v for v in values[p] if not np.isclose(v, base[p])], "osD^vPX")) for p in SWEEP_COLOR}
    rows = (len(names) + 1) // 2
    fig = plt.figure(figsize=(13, 3.2 * rows + 1.6), layout="constrained")
    panels, strip = fig.subfigures(2, 1, height_ratios=[3.2 * rows, 1.6])
    for letter, sf, s in zip("abcdefghij", panels.subfigures(rows, 2, wspace=0.15).flat, names):
        sf.suptitle(LABELS[s], fontweight="bold")
        sf.text(0.01, 0.98, f"({letter})", ha="left", va="top", transform=sf.transSubfigure, **PANEL_LABEL)
        sub = d[d["scheme"] == s]
        for ax, (col, label, title, pred) in zip(sf.subplots(1, len(QUANTITIES)), QUANTITIES):
            ax.set(xlabel=f"measured {label}", ylabel=f"predicted {label}")
            v = sub[[f"{col}_mean", f"{pred}_mean"]].dropna()
            if (v == 0).all().all():  # context-free schemes: Delta = 0 exactly, measured and predicted
                ax.text(0.5, 0.5, "0 by definition", transform=ax.transAxes, ha="center", va="center", color="#8a8985")
                ax.set(xticks=[], yticks=[])
                continue
            for _, r in sub.iterrows():
                p = r["sweep"]
                style = dict(marker="*", color=BASE_COLOR, ms=12, zorder=3) if p == "base" else \
                    dict(marker=markers[p][r[p]], color=SWEEP_COLOR[p], ms=8, zorder=2)
                ax.plot(r[f"{col}_mean"], r[f"{pred}_mean"], mec="white", mew=0.6, ls="none", **style)
            lo, hi = v.min().min(), v.max().max()
            if col == "H" or (col == "Delta" and lo > 0):
                lim = (lo / 1.2, hi * 1.2)
                log_axes(ax, np.log10(hi / lo))
            else:
                pad = 0.05 * (hi - lo) if hi > lo else 0.01
                lim = (lo - pad, hi + pad)
                linear_ticks(ax)
            ax.plot(lim, lim, color="#8a8985", lw=1, ls="--", zorder=0)  # y = x
            ax.set(xlim=lim, ylim=lim)
    # legends: one per sweep, titled with its name and aligned at the top, its values in order in columns of 3;
    # the base value is the star
    def entry(p, v):
        label = f"{v / 1e6:g} Mb" if p == "n" else f"{v:g}"
        if np.isclose(v, base[p]):
            return plt.Line2D([], [], ls="none", marker="*", color=BASE_COLOR, ms=12, label=label)
        return plt.Line2D([], [], ls="none", marker=markers[p][v], color=SWEEP_COLOR[p], ms=8, label=label)

    ncols = [-(-len(values[p]) // 3) for p in SWEEP_COLOR]
    for ax, p, ncol in zip(strip.subplots(1, len(SWEEP_COLOR), width_ratios=ncols), SWEEP_COLOR, ncols):
        ax.axis("off")
        legend = ax.legend(handles=[entry(p, v) for v in values[p]], title=f"{SWEEP_NAME[p]} sweep", ncol=ncol,
                           loc="upper center", columnspacing=1.0, handletextpad=0.2)
        legend.get_title().set_fontweight("bold")
    save(fig, f"{out}/schemes")


def scheme_table(summ, out):
    """Per scheme: measured and predicted values at the base setting, and how far the predictions are off over every
    setting (H: relative error; Delta: absolute error)."""
    rows = []
    for s in schemes_in(summ):
        sub = summ[summ["scheme"] == s]
        b = sub[sub["sweep"] == "base"].iloc[0]
        r = {"scheme": s, "settings": len(sub)}
        for c in ["d", "H", "H_pred", "H_true", "Delta", "Delta_pred", "Ehits", "U_N"]:
            r[c] = b[f"{c}_mean"]
        h_err = (sub["H_pred_mean"] / sub["H_mean"] - 1).abs()
        r["H_rel_err_median"], r["H_rel_err_max"] = h_err.median(), h_err.max()
        r["Delta_abs_err_median"] = (sub["Delta_pred_mean"] - sub["Delta_mean"]).abs().median()
        rows.append(r)
    pd.DataFrame(rows).to_csv(f"{out}/schemes.csv", index=False, float_format="%.4g")


# the table's schemes, in its order: minimizers, then context-free schemes, then the references
TABLE_ORDER = ["random-min", "lex-min", "weighted-min", "mod-min", "open-sync", "closed-sync", "fracminhash", "modmer",
               "all"]


def scheme_latex(summ, out):
    """schemes.tex: per scheme, measured and predicted H and Delta at the base setting (mean over genomes and hash
    seeds) and the relative error of the predicted H over all settings (median / largest); --- where a scheme has
    no prediction."""
    def number(x, digits):  # 4 decimals for H; 2 significant digits for Delta; 0 as 0
        if pd.isna(x):
            return "---"
        if x == 0:
            return "0"
        if digits == "H":
            return f"{x:.4f}"
        return np.format_float_positional(x, precision=2, unique=False, fractional=False, trim="-")

    def percent(x):
        return "---" if pd.isna(x) else f"{100 * x:.1f}\\%" if x < 0.1 else f"{100 * x:.0f}\\%"

    lines = []
    for s in [s for s in TABLE_ORDER if s in set(summ["scheme"])]:
        sub = summ[summ["scheme"] == s]
        b = sub[sub["sweep"] == "base"].iloc[0]
        err = (sub["H_pred_mean"] / sub["H_mean"] - 1).abs()
        error = "---" if err.isna().all() else f"{percent(err.median())} / {percent(err.max())}"
        lines.append(f"{LABELS[s].replace('k-mers', '$k$-mers')} & {number(b['H_mean'], 'H')} & {number(b['H_pred_mean'], 'H')} & "
                     f"{number(b['Delta_mean'], 'D')} & {number(b['Delta_pred_mean'], 'D')} & {error} \\\\")
    b = summ[summ["sweep"] == "base"].iloc[0]
    settings = summ.groupby("scheme").size().max()
    base = f"$n={b['n'] / 1e6:g}$\\,Mb, $\\alpha={b['alpha']:g}$, $\\theta={b['theta']:g}$, $k={b['k']:g}$, $w={b['w']:g}$"
    note = (f"Base setting: {base}; mean of {b['n_runs']:g} runs. Last column: relative error of the predicted "
            f"$\\Hu_\\sch$ over all {settings} settings (median / largest).")
    tex = rf"""\begin{{table}}[t]
\centering
\caption{{Measured and predicted uniqueness entropy $\Hu_\sch$ and false uniqueness $\Delta_\sch$ at the base setting, and the error of the predicted $\Hu_\sch$ across settings.}}
\label{{tab:sim-schemes}}
\small
\begin{{tabular}}{{lrrrrc}}
\toprule
 & \multicolumn{{2}}{{c}}{{$\Hu_\sch$}} & \multicolumn{{2}}{{c}}{{$\Delta_\sch$}} & error of predicted $\Hu_\sch$ \\
\cmidrule(lr){{2-3}} \cmidrule(lr){{4-5}}
scheme & measured & predicted & measured & predicted & median / largest \\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\par\smallskip
\begin{minipage}{\linewidth}\footnotesize
""" + note + r"""
\end{minipage}
\end{table}
"""
    with open(f"{out}/schemes.tex", "w") as fh:
        fh.write(tex)


# ---------------------------------------------------------------------------
# real genomes
# ---------------------------------------------------------------------------
def real():
    """ue's tables, summary.csv (mean and sd over hash seeds, with the predictions' relative errors) and
    predicted against measured H and Δ per genome and scheme."""
    panel = load("results/real/*.json", "panel")
    if panel.empty:
        return
    out = outdir("real")
    for t in ("spectrum_summary", "k_sweep", "alpha_from_sketch", "copies", "context_kernels"):
        load("results/real/*.json", t).to_csv(f"{out}/{t}.csv", index=False)
    panel.to_csv(f"{out}/raw.csv", index=False)
    summ = summarise(panel, ["genome"])
    for p in [*PRED, "H_pred_measured_d"]:
        summ[f"{p}_rel_err"] = summ[f"{p}_mean"] / summ["H_mean"] - 1
    summ.to_csv(f"{out}/summary.csv", index=False)

    with plt.rc_context(WIDE_FONTS):
        real_figure(summ, out)


def real_figure(summ, out):
    """Predicted against measured H (three predictions) and Delta, one row of panels (a)-(d), one point per genome and
    scheme (mean over hash seeds); color: scheme, marker: genome."""
    d = summ[summ["H_pred_mean"].notna()]
    genomes = [g for g in GENOMES if g in set(d["genome"])] + sorted(set(d["genome"]) - set(GENOMES))
    marker = dict(zip(genomes, "osD^vPXh"))
    panels = [(p, "H", title) for p, title in PRED.items()] + [("Delta_pred", "Delta", "false uniqueness Δ, plug-in")]
    fig = plt.figure(figsize=(17, 6.2), layout="constrained")
    plots, strip = fig.subfigures(2, 1, height_ratios=[4.9, 1.3])
    for letter, sf, (pcol, mcol, title) in zip("abcd", plots.subfigures(1, 4, wspace=0.06), panels):
        sf.suptitle(f"({letter})", x=0.01, ha="left", **PANEL_LABEL)
        ax = sf.subplots()
        for _, r in d.iterrows():
            ax.plot(r[f"{mcol}_mean"], r[f"{pcol}_mean"], marker[r["genome"]], color=COLORS[r["scheme"]], ms=8,
                    mec="white", mew=0.5)
        v = d[[f"{mcol}_mean", f"{pcol}_mean"]]
        lo, hi = v.min().min(), v.max().max()
        if mcol == "H":  # H spans bacteria to plants
            lim = (lo / 1.2, hi * 1.2)
            log_axes(ax, np.log10(hi / lo))
        else:
            pad = 0.05 * (hi - lo) if hi > lo else 0.01
            lim = (lo - pad, hi + pad)
            linear_ticks(ax)
        ax.plot(lim, lim, color="#8a8985", lw=1, ls="--", zorder=0)  # y = x
        q = "H" if mcol == "H" else "Δ"
        ax.set(xlim=lim, ylim=lim, xlabel=f"measured {q} (bits)", ylabel=f"predicted {q} (bits)", title=title)
    # legends under the panels: the genomes in one line, then the schemes
    ax = strip.subplots()
    ax.axis("off")
    genome_line = [plt.Line2D([], [], ls="none", label="genome")]
    genome_line += [plt.Line2D([], [], color="#0b0b0b", marker=marker[g], mfc="none", ls="none", ms=8,
                               label=GENOMES.get(g, g)) for g in genomes]
    first = ax.legend(handles=genome_line, loc="upper center", ncol=len(genome_line), columnspacing=1.0,
                      handletextpad=0.3)
    first.get_texts()[0].set_fontweight("bold")
    ax.add_artist(first)
    schemes = [plt.Line2D([], [], color=COLORS[s], marker="o", ls="none", ms=8, label=LABELS[s]) for s in schemes_in(d)]
    ax.legend(handles=schemes, loc="lower center", ncol=len(schemes), columnspacing=1.2, handletextpad=0.3)
    save(fig, f"{out}/predicted_vs_measured")


def main():
    os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))  # the repository
    simulated()
    real()


if __name__ == "__main__":
    main()
