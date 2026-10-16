# merged_greedy_plotting.py
#
# Single merged figure for Supplemental Fig. "greedy_baseline" (Reviewer 3,
# comment 7): one shared y axis, one shared x axis (via supxlabel), and the
# missing-rate legend embedded in the empty band of panel A, so the axis
# labels are drawn once instead of once per panel and can be set larger.
#
#   A. CoEvoLink, supervised setting, with genomic information
#      (reproduces figures_boxplots_final/normal_f1_genomic.pdf)
#   B. Degree-product greedy baseline (virus degree x host degree),
#      same instances, same budget -- F1 rescaled by a constant factor
#      (baseline_plotting.PROD_F1_CEILING) for visual comparability across
#      missing rates; see baseline_plotting.py for the factor and the raw,
#      unscaled numbers.
#
# The figure is authored at its final printed size (6.5 in = \textwidth of
# the 11pt article class with 1in margins) and included at width=\textwidth,
# so it renders at 1:1 and the point sizes below are literally what appears
# on paper.
#
# MATCH_INSTANCES=1 additionally restricts both panels to the configurations
# they have in common. That is off by default because it is not a neutral
# change: it drops the 50%-missing runs where CoEvoLink scores zero, which
# tend to be where the greedy baseline scores lowest too, so panel B's 50%
# boxes move up noticeably under it.
#
# Panels share the y axis, so the y label and its tick labels are drawn once;
# the x label is shared via supxlabel and the legend sits in the empty band of
# panel A.  That recovered space is what pays for the larger type.

import importlib.util
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# ==============================
# Configuration
# ==============================
FIGDIR = "figures_merged"
STEM = "greedy_prod_baseline_merged"

# restrict both panels to their common configurations (see header)
MATCH_INSTANCES = os.environ.get("MATCH_INSTANCES", "0") == "1"

# jittered per-replicate points: legible at 10x10 in, muddy at 6.5x3.2
SHOW_POINTS = os.environ.get("SHOW_POINTS", "0") == "1"

FLIP_MAP = {1000: "10%", 2000: "20%", 10000: "50%"}
FLIP_LEVELS = ["10%", "20%", "50%"]
PALETTE = ["#aec7e8", "#1f77b4", "#08306b"]  # light -> dark blue
MAX_RATE = 0.5
KEYS = ["seed", "corrupt", "r01", "r10"]

# Point sizes as printed, since the figure is drawn at final size.
# Bumped up from the first cut (9.5/11.0/11.0/9.0): sharing one y label and
# one x label between the two panels frees the width each used to spend on
# its own copy, so the type can grow and still fit at 1:1 on the page.
FS_TICK = 11.0
FS_LABEL = 13.5
FS_PANEL = 13.5
FS_LEGEND = 10.5

# reuse the loaders in baseline_plotting.py rather than duplicating them
_spec = importlib.util.spec_from_file_location("baseline_plotting",
                                               "baseline_plotting.py")
bp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bp)


# ==============================
# Data
# ==============================
def subset(df):
    """The same slice of the grid the existing figures use."""
    df = df[(df["r01"] == df["r10"]) & (df["r01"] <= MAX_RATE)].copy()
    df["rate"] = df["r01"]
    df = df[df["corrupt"].isin(FLIP_MAP)].copy()
    df["flips"] = df["corrupt"].map(FLIP_MAP)
    return df


def load():
    """CoEvoLink and the degree-product greedy baseline, same configurations."""
    # same source and same filter as normal_f1_genomic.pdf (box_plotting.py)
    coevolink = pd.read_csv("simulation_all_metrics_genomic.csv")
    coevolink = coevolink[~((coevolink["precision"] <= 0.0)
                            | (coevolink["recall"] <= 0.0)
                            | (coevolink["f1"] <= 0.0))].copy()
    greedy = bp.collect(bp.BASELINE_ROOT, "supervised_greedy_prod_metrics.json")

    coevolink, greedy = subset(coevolink), subset(greedy)
    if MATCH_INSTANCES:
        common = coevolink[KEYS].merge(greedy[KEYS], on=KEYS, how="inner")
        coevolink = coevolink.merge(common, on=KEYS, how="inner")
        greedy = greedy.merge(common, on=KEYS, how="inner")

    # same rescaling as greedy_prod_normal_f1.pdf (baseline_plotting.py):
    # multiply F1 by a constant factor so the highest value maps to
    # PROD_F1_CEILING, for visual comparability across missing rates.
    prod_scale = bp.PROD_F1_CEILING / greedy["f1"].max()
    greedy["f1_raw"] = greedy["f1"]
    greedy["f1"] = greedy["f1"] * prod_scale

    print(f"match_instances={MATCH_INSTANCES}  "
          f"rows: CoEvoLink {len(coevolink)}, greedy_prod {len(greedy)}")
    print(f"greedy_prod scaled by {prod_scale:.4f} "
          f"(max raw F1 {greedy['f1_raw'].max():.4f} -> {bp.PROD_F1_CEILING})")
    print(f"Median F1  CoEvoLink {coevolink['f1'].median():.4f}"
          f"  greedy_prod (scaled) {greedy['f1'].median():.4f}"
          f"  greedy_prod (raw) {greedy['f1_raw'].median():.4f}")
    print("\nReplicates per (rate, missing rate):")
    print(coevolink.groupby(["rate", "flips"]).size().unstack()[FLIP_LEVELS])
    return coevolink, greedy


# ==============================
# Figure
# ==============================
def panel(ax, df, letter):
    sns.boxplot(
        x="rate", y="f1", hue="flips", data=df, hue_order=FLIP_LEVELS,
        palette=PALETTE, showfliers=False, ax=ax, width=0.62,
        linewidth=0.8, legend=False,
    )
    if SHOW_POINTS:
        sns.stripplot(
            x="rate", y="f1", hue="flips", data=df, hue_order=FLIP_LEVELS,
            dodge=True, alpha=0.35, size=2, linewidth=0,
            palette=PALETTE, ax=ax, legend=False,
        )

    ax.yaxis.grid(True, linestyle="--", linewidth=0.6, alpha=0.6)
    ax.xaxis.grid(False)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(labelsize=FS_TICK)
    # above the axes, not inside them: panel A's 10% box reaches ~0.97 and
    # collided with an in-axes letter
    ax.text(0.0, 1.015, letter, transform=ax.transAxes,
            fontsize=FS_PANEL, fontweight="bold", va="bottom", ha="left")


def main():
    coevolink, greedy = load()

    # baseline_plotting.py sets a 10x10-inch theme at import; start clean so
    # the point sizes above are the only thing driving type size here.
    plt.rcParams.update(plt.rcParamsDefault)
    sns.set_context("paper", font_scale=1.0)
    sns.set_style("ticks")

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(6.5, 3.2), sharey=True,
        gridspec_kw={"wspace": 0.06},
    )

    panel(ax_a, coevolink, "A")
    panel(ax_b, greedy, "B")

    ax_a.set_ylabel("F1 Score", fontsize=FS_LABEL)
    fig.supxlabel(r"Transition rate $r$", fontsize=FS_LABEL, y=0.02)

    # legend in the empty band of panel A, which no CoEvoLink box reaches
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="0.25",
                             linewidth=0.8) for c in PALETTE]
    ax_a.legend(
        handles, FLIP_LEVELS, title="Missing rate",
        loc="lower left", bbox_to_anchor=(0.01, 0.01), ncol=3,
        fontsize=FS_LEGEND, title_fontsize=FS_LEGEND, frameon=False,
        handlelength=1.1, handletextpad=0.5, columnspacing=1.0,
        borderpad=0.2, labelspacing=0.3,
    )

    sns.despine(fig=fig)
    fig.patch.set_facecolor("white")
    # not tight_layout: it does not understand supxlabel and warns; the
    # margins below plus bbox_inches="tight" on save do the same job.
    fig.subplots_adjust(left=0.085, right=0.995, top=0.92, bottom=0.17)

    os.makedirs(FIGDIR, exist_ok=True)
    suffix = ("_points" if SHOW_POINTS else "")
    suffix += ("_matched" if MATCH_INSTANCES else "")
    path = os.path.join(FIGDIR, STEM + suffix)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{path}.{ext}", bbox_inches="tight",
                    dpi=300 if ext != "svg" else None)
    plt.close(fig)
    print(f"\nSaved: {path}.{{svg,pdf,png}}")


if __name__ == "__main__":
    main()
