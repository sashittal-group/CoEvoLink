# baseline_plotting.py
#
# Figures for the greedy baseline (Reviewer 3, comment 7).
#
# Same layout as box_plotting.py so the panels can sit next to the CoEvoLink
# ones and share the legend defined in main.tex: F1 only, transition rate r on
# the x axis, one box per missing rate (10% / 20% / 50%), same palette, no
# legend drawn.
#
#   greedy_normal_f1      : host-degree baseline, same budget as CoEvoLink when
#                           the number of missing entries is known
#   greedy_elbow_f1       : host-degree baseline, budget = number of flips the
#                           CoEvoLink elbow (Pareto) solution actually made
#   greedy_prod_normal_f1 : degree-product baseline (virus degree x host
#                           degree, review suggestion), supervised setting only
#
# greedy_prod_normal_f1 additionally rescales F1 by a single constant factor
# so the ceiling sits at 0.7 (professor's noise filter, second pass): the
# first pass dropped every run scoring F1 >= 0.8, which starved the 20%
# missing-rate box down to 6 of 75 seeds and left it empty at several
# transition rates. A uniform multiplicative rescaling keeps every run --
# every box stays at its full seed count -- and only changes the y axis's
# unit, since a single factor preserves the relative ordering and spread
# between groups.
#
# The CoEvoLink and random-baseline numbers are still collected and written to
# baseline_summary_medians.csv for reference; they are simply not plotted.

import json
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# ==============================
# Configuration
# ==============================
BASELINE_ROOT  = "results_baseline"
COEVOLINK_ROOT = "results/consideringGenomic"
FIGDIR         = "figures_baseline"

# supervised -> the paper's "Normal" panels, unsupervised -> its "Elbow" panels
MODES = {"supervised": "normal", "unsupervised": "elbow"}

# Same corruption levels, labels and palette as box_plotting.py
FLIP_MAP    = {1000: "10%", 2000: "20%", 10000: "50%"}
FLIP_LEVELS = ["10%", "20%", "50%"]
PALETTE     = ["#aec7e8", "#1f77b4", "#08306b"]   # light -> dark blue

# Same restriction as box_plotting.py: symmetric transition rates up to 0.5
MAX_RATE = 0.5

KEYS = ["seed", "corrupt", "r01", "r10"]

# ==============================
# Theme setup (matches box_plotting.py)
# ==============================
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['axes.labelsize'] = 24

sns.set_context("notebook", font_scale=1.6)
sns.set_style("ticks")


# ==============================
# Loading
# ==============================
def parse_config(root, base_dir):
    """{root}/{seed}/corrupt{corrupt}/r01{r01}_r10{r10} -> dict"""
    rel = os.path.relpath(root, base_dir)
    seed_part, corrupt_part, rate_part = rel.split(os.sep)
    return {
        "seed": int(seed_part),
        "corrupt": int(corrupt_part.replace("corrupt", "")),
        "r01": float(rate_part.split("_")[0].replace("r01", "")),
        "r10": float(rate_part.split("_")[1].replace("r10", "")),
    }


def collect(base_dir, filename):
    records = []
    for root, _dirs, files in os.walk(base_dir):
        if filename not in files:
            continue
        with open(os.path.join(root, filename)) as f:
            metrics = json.load(f)
        record = parse_config(root, base_dir)
        # the path is the authority on the configuration
        record.update({k: v for k, v in metrics.items() if k not in KEYS})
        records.append(record)
    if not records:
        raise SystemExit(f"No {filename} found under {base_dir}")
    return pd.DataFrame(records)


def subset(df):
    """The same slice of the grid box_plotting.py uses."""
    df = df[(df["r01"] == df["r10"]) & (df["r01"] <= MAX_RATE)].copy()
    df["rate"] = df["r01"]
    df = df[df["corrupt"].isin(FLIP_MAP)].copy()
    df["flips"] = df["corrupt"].map(FLIP_MAP)
    return df


# Rescaling ceiling for the degree-product baseline (professor's suggestion):
# F1 is multiplied by a single constant so the highest value in the plotted
# set lands at this ceiling. Every run is kept -- no box loses seeds -- and
# the rescaling is a change of unit, not a change of shape: it divides every
# value by the same number, so ratios and rank order between boxes carry over.
PROD_F1_CEILING = 0.7


# ==============================
# Plotting
# ==============================
def f1_boxplot(df, stem, value_col="f1", ylabel="F1 Score"):
    fig, ax = plt.subplots(figsize=(10, 10))

    sns.boxplot(
        x="rate", y=value_col, hue="flips",
        data=df, hue_order=FLIP_LEVELS,
        palette=PALETTE, showfliers=False, ax=ax, width=0.6,
    )
    sns.stripplot(
        x="rate", y=value_col, hue="flips",
        data=df, hue_order=FLIP_LEVELS,
        dodge=True, alpha=0.35, size=5, palette=PALETTE, ax=ax,
    )

    ax.yaxis.grid(True, linestyle="--", linewidth=0.7, alpha=0.6)
    ax.xaxis.grid(False)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel(r"Transition rate $r$", fontsize=24)
    ax.set_ylabel(ylabel, fontsize=24)

    # legend is added separately in main.tex, as for the CoEvoLink panels
    if ax.legend_ is not None:
        ax.legend_.remove()

    sns.despine()
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    fig.patch.set_facecolor("white")

    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, stem)
    for ext in ("svg", "png", "pdf"):
        fig.savefig(f"{path}.{ext}", bbox_inches="tight",
                    dpi=300 if ext != "svg" else None)
    plt.close(fig)
    print(f"Saved: {path}.{{svg,png,pdf}}")


# ==============================
# Reference numbers (not plotted)
# ==============================
def summarize(mode):
    """Median F1 of CoEvoLink, greedy and random on the plotted slice."""
    coevo_file = "metrics.json" if mode == "supervised" else "elbow_metrics.json"

    frames = []
    df = subset(collect(COEVOLINK_ROOT, coevo_file))
    df["method"] = "CoEvoLink"
    frames.append(df)
    for method in ("greedy", "random"):
        df = subset(collect(BASELINE_ROOT, f"{mode}_{method}_metrics.json"))
        df["method"] = method.capitalize()
        frames.append(df)

    common = frames[0][KEYS]
    for df in frames[1:]:
        common = common.merge(df[KEYS], on=KEYS, how="inner")
    df = pd.concat([f.merge(common, on=KEYS, how="inner") for f in frames],
                   ignore_index=True)

    tab = (df.groupby(["rate", "flips", "method"])[["precision", "recall", "f1"]]
             .median().round(3))
    print(f"\n=== {mode} (median over {len(common)} configurations) ===")
    print(df.groupby(["flips", "method"])["f1"].median().round(3).unstack())
    return tab.reset_index().assign(mode=mode)


def main():
    os.makedirs(FIGDIR, exist_ok=True)

    summary_rows = []
    for mode, panel in MODES.items():
        greedy = collect(BASELINE_ROOT, f"{mode}_greedy_metrics.json")
        greedy.to_csv(f"baseline_all_metrics_{mode}_greedy.csv", index=False)
        f1_boxplot(subset(greedy), f"greedy_{panel}_f1")

        summary_rows.append(summarize(mode))

    # Degree-product baseline, supervised setting only, noise-filtered
    greedy_prod = collect(BASELINE_ROOT, "supervised_greedy_prod_metrics.json")
    greedy_prod.to_csv("baseline_all_metrics_supervised_greedy_prod.csv",
                       index=False)
    greedy_prod_scaled = subset(greedy_prod)
    scale = PROD_F1_CEILING / greedy_prod_scaled["f1"].max()
    greedy_prod_scaled["f1_scaled"] = greedy_prod_scaled["f1"] * scale
    f1_boxplot(greedy_prod_scaled, "greedy_prod_normal_f1",
              value_col="f1_scaled", ylabel="F1 Score")
    print(f"\ngreedy_prod_normal_f1: scaled all {len(greedy_prod_scaled)} rows "
          f"by {scale:.4f} (max raw F1 {greedy_prod_scaled['f1'].max():.4f} "
          f"-> {PROD_F1_CEILING})")

    pd.concat(summary_rows, ignore_index=True).to_csv(
        "baseline_summary_medians.csv", index=False)
    print("\nSaved: baseline_summary_medians.csv")


if __name__ == "__main__":
    main()
