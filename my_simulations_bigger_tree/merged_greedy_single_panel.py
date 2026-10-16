# merged_greedy_single_panel.py
#
# The one-axes alternative to merged_greedy_plotting.py: instead of two
# panels, both methods share a single set of axes, so every transition rate r
# holds two sub-blocks of three boxes -- CoEvoLink then the greedy baseline,
# one box per missing rate.  Thirty boxes in all.
#
# The two sub-blocks are separated by alternating background bands (CoEvoLink
# on white, the baseline on a light tint) and named on the inner row of the x
# axis, so the pairing is readable without decoding the colours.  Colour still
# carries the missing rate alone, exactly as in the two-panel version.
#
# Data handling is imported from merged_greedy_plotting, so the sources and
# the zero-metric filter of the original scripts are unchanged: each method
# contributes exactly the rows its own figure plots.

import importlib.util
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.transforms import blended_transform_factory

_spec = importlib.util.spec_from_file_location("merged_greedy_plotting",
                                               "merged_greedy_plotting.py")
mg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mg)

FIGDIR = "figures_merged"
STEM = "greedy_baseline_single"

FLIP_LEVELS = mg.FLIP_LEVELS
PALETTE = mg.PALETTE

METHODS = ["CoEvoLink", "Greedy"]
BAND_COLORS = {"CoEvoLink": "none", "Greedy": "#eef1f5"}

# 30 boxes is dense enough that the per-replicate dots turn to mush
SHOW_POINTS = os.environ.get("SHOW_POINTS", "0") == "1"

FS_TICK = 9.5      # transition rates, y ticks
FS_INNER = 8.5     # method names under each sub-block
FS_LABEL = 11.0
FS_LEGEND = 9.0


def build():
    coevolink, greedy = mg.load()
    coevolink = coevolink.assign(method="CoEvoLink")
    greedy = greedy.assign(method="Greedy")
    df = pd.concat([coevolink, greedy], ignore_index=True)
    df["slot"] = df["rate"].astype(str) + "|" + df["method"]
    rates = sorted(df["rate"].unique())
    order = [f"{r}|{m}" for r in rates for m in METHODS]
    return df, rates, order


def main():
    df, rates, order = build()

    plt.rcParams.update(plt.rcParamsDefault)
    sns.set_context("paper", font_scale=1.0)
    sns.set_style("ticks")

    fig, ax = plt.subplots(figsize=(6.5, 4.0))

    # background bands first, so boxes and grid sit on top of them
    for i, slot in enumerate(order):
        color = BAND_COLORS[slot.split("|")[1]]
        if color != "none":
            ax.axvspan(i - 0.5, i + 0.5, color=color, lw=0, zorder=0)
    # thin rule between transition rates
    for i in range(2, len(order), 2):
        ax.axvline(i - 0.5, color="0.75", lw=0.6, zorder=0.5)

    sns.boxplot(
        x="slot", y="f1", hue="flips", data=df, order=order,
        hue_order=FLIP_LEVELS, palette=PALETTE, showfliers=False,
        ax=ax, width=0.78, linewidth=0.7, legend=False, zorder=3,
    )
    if SHOW_POINTS:
        sns.stripplot(
            x="slot", y="f1", hue="flips", data=df, order=order,
            hue_order=FLIP_LEVELS, dodge=True, alpha=0.35, size=1.6,
            linewidth=0, palette=PALETTE, ax=ax, legend=False, zorder=4,
        )

    ax.yaxis.grid(True, linestyle="--", linewidth=0.6, alpha=0.6, zorder=1)
    ax.xaxis.grid(False)
    ax.set_ylim(0, 1.02)
    ax.set_xlim(-0.5, len(order) - 0.5)
    ax.set_ylabel("F1 Score", fontsize=FS_LABEL)
    ax.tick_params(axis="y", labelsize=FS_TICK)

    # inner x row: the method under each sub-block
    ax.set_xlabel("")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([s.split("|")[1] for s in order], fontsize=FS_INNER)
    ax.tick_params(axis="x", length=0, pad=2)

    # outer x row: the transition rate, centred under each pair.  Offsets are
    # in points, not axes fractions, so they hold if the figure is resized.
    tform = blended_transform_factory(ax.transData, ax.transAxes)
    for i, r in enumerate(rates):
        ax.annotate(f"{r}", xy=(2 * i + 0.5, 0), xycoords=tform,
                    xytext=(0, -(FS_INNER + 9)), textcoords="offset points",
                    ha="center", va="top", fontsize=FS_TICK)
    ax.annotate(r"Transition rate $r$", xy=(0.5, 0), xycoords=ax.transAxes,
                xytext=(0, -(FS_INNER + FS_TICK + 18)),
                textcoords="offset points",
                ha="center", va="top", fontsize=FS_LABEL)

    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="0.25",
                             linewidth=0.7) for c in PALETTE]
    ax.legend(
        handles, FLIP_LEVELS, title="Missing rate",
        loc="lower center", bbox_to_anchor=(0.5, 1.005), ncol=3,
        fontsize=FS_LEGEND, title_fontsize=FS_LEGEND, frameon=False,
        handlelength=1.1, handletextpad=0.5, columnspacing=1.2,
        borderpad=0.2,
    )

    sns.despine(fig=fig)
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.085, right=0.995, top=0.915, bottom=0.175)

    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, STEM + ("_points" if SHOW_POINTS else ""))
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{path}.{ext}", bbox_inches="tight",
                    dpi=300 if ext != "svg" else None)
    plt.close(fig)
    print(f"\nSaved: {path}.{{svg,pdf,png}}")


if __name__ == "__main__":
    main()
