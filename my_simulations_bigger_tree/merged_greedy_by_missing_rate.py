# merged_greedy_by_missing_rate.py
#
# The by-missing-rate alternative to merged_greedy_plotting.py: one panel per
# missing rate (A 10%, B 20%, C 50%) instead of one panel per method.  Within
# each panel every transition rate r holds two boxes side by side -- CoEvoLink
# and the degree-product greedy baseline -- so colour carries the method and
# the missing rate is written inside the panel.
#
# Laid out like main-text Figure 6: boxed panels in one row, the y label and
# the method legend drawn once (below the x label), panel letters above the
# top-left corners.
#
# Data handling is imported from merged_greedy_plotting, so the sources, the
# zero-metric filter, and the greedy F1 rescaling are unchanged.
#
# The figure is authored at its printed size (6.5 in = \textwidth of the 11pt
# article class with 1in margins) and included at width=\textwidth, so the
# point sizes below are literally what appears on paper.

import importlib.util
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

_spec = importlib.util.spec_from_file_location("merged_greedy_plotting",
                                               "merged_greedy_plotting.py")
mg = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mg)

FIGDIR = "figures_merged"
STEM = "greedy_prod_baseline_by_missing"

METHODS = ["CoEvoLink", "Greedy heuristic"]
METHOD_COLORS = ["#1f4e8c", "#a6a6a6"]  # blue / gray: separated by lightness

# jittered per-replicate points: off by default, as in the merged figure
SHOW_POINTS = os.environ.get("SHOW_POINTS", "0") == "1"

FIGSIZE = (6.5, 2.7)

# point sizes as printed, since the figure is drawn at final size
FS_TICK = 8.0
FS_LABEL = 9.5
FS_PANEL = 10.0
FS_LEGEND = 8.5
FS_INSET = 8.5


def build():
    coevolink, greedy = mg.load()
    df = pd.concat([coevolink.assign(method=METHODS[0]),
                    greedy.assign(method=METHODS[1])], ignore_index=True)
    return df


def panel(ax, df, level, letter):
    sns.boxplot(
        x="rate", y="f1", hue="method", data=df, hue_order=METHODS,
        palette=METHOD_COLORS, showfliers=False, ax=ax, width=0.62,
        linewidth=0.6, saturation=1, legend=False,
    )
    if SHOW_POINTS:
        sns.stripplot(
            x="rate", y="f1", hue="method", data=df, hue_order=METHODS,
            dodge=True, alpha=0.35, size=2, linewidth=0,
            palette=METHOD_COLORS, ax=ax, legend=False,
        )

    ax.yaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    ax.xaxis.grid(False)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(labelsize=FS_TICK, length=2.5, width=0.6, pad=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)

    # bottom right is empty in every panel: the lowest whisker (greedy at
    # 50%, r = 0.1) stops near 0.16 and sits at the left end
    ax.text(0.97, 0.03, f"Missing rate {level}", transform=ax.transAxes,
            fontsize=FS_INSET, ha="right", va="bottom")
    ax.text(0.0, 1.02, letter, transform=ax.transAxes,
            fontsize=FS_PANEL, fontweight="bold", ha="left", va="bottom")


def main():
    df = build()

    plt.rcParams.update(plt.rcParamsDefault)
    sns.set_context("paper", font_scale=1.0)
    sns.set_style("ticks")

    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE, sharey=True,
                             gridspec_kw={"wspace": 0.06})
    for ax, level, letter in zip(axes, mg.FLIP_LEVELS, "ABC"):
        panel(ax, df[df["flips"] == level], level, letter)

    axes[0].set_ylabel("F1 Score", fontsize=FS_LABEL, labelpad=3)
    fig.supxlabel(r"Transition rate $r$", fontsize=FS_LABEL, y=0.095)

    # one legend for the figure, below the x label: the method colours are
    # the same in all panels
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="0.25",
                             linewidth=0.6) for c in METHOD_COLORS]
    fig.legend(
        handles, METHODS, loc="lower center", bbox_to_anchor=(0.5, 0.0),
        ncol=2, fontsize=FS_LEGEND, frameon=False,
        handlelength=1.0, handleheight=0.7, handletextpad=0.4,
        columnspacing=1.5, borderpad=0.1,
    )

    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.07, right=0.995, top=0.92, bottom=0.215)

    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, STEM + ("_points" if SHOW_POINTS else ""))
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{path}.{ext}", bbox_inches="tight",
                    dpi=300 if ext != "svg" else None)
    plt.close(fig)
    print(f"Saved: {path}.{{svg,pdf,png}}")


if __name__ == "__main__":
    main()
