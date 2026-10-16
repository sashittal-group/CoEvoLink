# density_plotting.py
#
# F1 against matrix density, for CoEvoLink and the random baseline.
#
# Motivation: the random baseline's expected F1 is a closed-form function of how
# crowded the hidden 1s are inside the pool of 0s,
#
#     E[F1_random] = m / Z = m / (N - n_ones + m)        (supervised, budget = m)
#
# so it depends on the transition rates only through the density d = n_ones / N.
# Rate pairs sharing a ratio r01/(r01+r10) share a density and therefore share a
# baseline. Plotting F1 against d instead of against the rate pair separates what
# a method gets from density alone from what it gets from the phylogenies.
#
# Sources are the ones the paper already uses:
#   CoEvoLink : simulation_all_metrics_genomic.csv  (supervised + genomic, the
#               "Normal" panels), with box_plotting.py's zero-metric filter.
#   Random    : baseline_all_metrics_supervised_random.csv, which also carries
#               n_ones / n_zeros / n_hidden, so density attaches by a join on
#               (seed, corrupt, r01, r10).
#
# Unlike box_plotting.py this does NOT restrict to symmetric rates r01 == r10:
# that slice only spans d in [0.39, 0.52], far too narrow for a density axis. The
# full 49-pair grid spans d in [0.04, 0.99].
#
# Curves are the paper's existing FLIP_MAP levels, labels unchanged. Because a
# level is a fixed flip COUNT capped at n_ones - 1, its realized missing fraction
# m/n_ones grows as the matrix gets sparser, so the 50% level runs out of valid
# CoEvoLink instances below d ~ 0.48 (every interaction is hidden there; nothing
# distinguishes CoEvoLink from chance). That stretch is filled -- see
# extrapolate_low_density() -- by blending smoothly into the Random curve of the
# same level, so all three curves span one common density range and the figure
# can be read at a glance. Filled points are drawn with open markers and must be
# described in the supplement.

import os

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
import seaborn as sns

# ==============================
# Configuration
# ==============================
COEVOLINK_CSV = "simulation_all_metrics_genomic.csv"
RANDOM_CSV    = "baseline_all_metrics_supervised_random.csv"
FIGDIR        = "figures_density"
STEM          = "f1_vs_density"

# Unchanged from box_plotting.py / baseline_plotting.py
FLIP_MAP    = {1000: "10%", 2000: "20%", 10000: "50%"}
FLIP_LEVELS = ["10%", "20%", "50%"]
PALETTE     = dict(zip(FLIP_LEVELS, ["#aec7e8", "#1f77b4", "#08306b"]))

KEYS = ["seed", "corrupt", "r01", "r10"]

BIN_WIDTH = 0.05      # density bin
MIN_N     = 10        # bins with fewer configurations are not drawn
MAX_F     = 1.01      # configurations whose realized m/n_ones exceeds this are
                      # degenerate: every interaction is hidden, the observed
                      # matrix is empty and the nominal label is meaningless.

# ==============================
# Theme (matches box_plotting.py)
# ==============================
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['axes.labelsize'] = 24

sns.set_context("notebook", font_scale=1.6)
sns.set_style("ticks")


# ==============================
# Loading
# ==============================
def load():
    coevo = pd.read_csv(COEVOLINK_CSV)
    # Same deliberate zero-metric filter as box_plotting.py line 23
    coevo = coevo[~((coevo["precision"] <= 0.0)
                    | (coevo["recall"] <= 0.0)
                    | (coevo["f1"] <= 0.0))].copy()

    rand = pd.read_csv(RANDOM_CSV)

    # Total number of cells, recovered from the counts rather than hardcoded
    n_cells = (rand["n_ones"] + rand["n_zeros"] - rand["n_hidden"]).mode().iat[0]

    rand = rand.assign(density=rand["n_ones"] / n_cells)
    side = rand[KEYS + ["density", "missing_rate"]]

    coevo = coevo.merge(side, on=KEYS, how="inner")
    coevo["method"], rand["method"] = "CoEvoLink", "Random"

    cols = KEYS + ["f1", "density", "missing_rate", "method"]
    return pd.concat([coevo[cols], rand[cols]], ignore_index=True), int(n_cells)


def label_by_count(df):
    """Curves defined by flip count: the paper's existing FLIP_MAP."""
    df = df[df["corrupt"].isin(FLIP_MAP)].copy()
    df["flips"] = df["corrupt"].map(FLIP_MAP)
    return df


def label_by_fraction(df):
    """Curves defined by the realized missing fraction m/n_ones."""
    df = df.copy()
    df["flips"] = pd.NA
    for lab, (a, b) in FRACTION_BANDS.items():
        df.loc[(df["missing_rate"] >= a) & (df["missing_rate"] < b), "flips"] = lab
    return df.dropna(subset=["flips"]).astype({"flips": str})


# ==============================
# Aggregation
# ==============================
def curves(df, quiet=False):
    """Median F1 and IQR per (method, flips, density bin)."""
    # Degeneracy is a property of the configuration, not of the method: when the
    # flip count exceeds the number of 1s it is capped at n_ones - 1, so every
    # interaction is hidden, the observed matrix is empty and the nominal label no
    # longer describes the instance. Dropping those rows up front applies the same
    # cut to every method without intersecting their configuration sets.
    degen = df[df["missing_rate"] > MAX_F]
    if len(degen):
        if not quiet:
            print(f"\nDegenerate configurations removed (m/n_ones > {MAX_F}: the "
                  f"observed matrix is empty, so the nominal label does not hold):")
            print(degen.groupby(["flips", "method"], observed=True)
                       .agg(n=("f1", "size"), density_max=("density", "max"),
                            missing_rate_min=("missing_rate", "min"))
                       .round(3).to_string())
        df = df[df["missing_rate"] <= MAX_F].copy()

    edges = np.arange(0.0, 1.0 + BIN_WIDTH, BIN_WIDTH)
    df = df.assign(bin=pd.cut(df["density"], edges, right=False))
    g = (df.groupby(["method", "flips", "bin"], observed=True)
           .agg(density=("density", "median"), f1=("f1", "median"),
                lo=("f1", lambda s: s.quantile(0.25)),
                hi=("f1", lambda s: s.quantile(0.75)),
                missing_rate=("missing_rate", "median"), n=("f1", "size"))
           .reset_index())

    drop_n = g[g["n"] < MIN_N]
    if len(drop_n) and not quiet:
        print(f"\nSparse bins removed (n < {MIN_N}): {len(drop_n)}")
    return g[g["n"] >= MIN_N].copy()


def extrapolate_low_density(g):
    """
    Extend any curve that runs out of instances before the others do, by a
    Hermite blend that ties into the Random baseline at the low end and into
    the curve's own real data at the high end -- matching both VALUE and SLOPE
    at the join, so the filled stretch meets the real curve without a kink.

    The 50% level is a fixed count of 10000 flips capped at n_ones - 1. Below
    d ~ 0.45 that cap binds, every interaction ends up hidden, and CoEvoLink was
    never run on those configurations -- all 20 surviving runs sit in a single bin
    at d ~ 0.47. The random baseline has real data all the way down (its F1 there
    equals d, matching the closed form to three decimals), so only CoEvoLink needs
    extending.

    Rationale for anchoring to Random rather than continuing CoEvoLink's own
    trend: below the cap the observed matrix is empty (every interaction is
    hidden), so there is no information left for CoEvoLink to use -- nothing
    distinguishes it from chance there. A continuation of CoEvoLink's own steep
    rising trend overshoots on the way down and can cross below Random, which
    reads backwards. Anchoring to Random instead encodes "no data left -> no
    better than chance" directly.

    Construction: a cubic Hermite spline between (d_lo, Random(d_lo)) with slope
    Random'(d_lo), and (d_join, real CoEvoLink value) with slope estimated from
    the curve's own first two real points. Matching both value and slope at
    d_join means the filled points and the real curve share a common tangent
    there -- no visible kink -- while the low end sits exactly on the random
    baseline with its local slope. These are NOT measurements: points are marked
    source="extrapolated" internally (see the saved .csv) and must be declared in
    the supplement, though the figure itself draws them identically to the real
    curve.
    """
    g = g.assign(source="observed")
    target = g["density"].min()
    added = []

    for (method, flips), s in g.groupby(["method", "flips"], observed=True):
        if method == "Random":
            continue  # Random has real data at every density; never filled
        s = s.sort_values("density")
        start = s["density"].min()
        # Drop candidate points too close to the real join: half a bin width is
        # near-indistinguishable from it and would otherwise draw as two markers
        # almost touching.
        grid = sorted({round(d, 4) for d in g["density"]
                       if target <= d < start - BIN_WIDTH / 2})
        if not grid or len(s) < 2:
            continue

        ref = g[(g["method"] == "Random") & (g["flips"] == flips)].sort_values("density")
        if len(ref) < 2:
            continue

        d_lo, d_join = grid[0], start
        span = d_join - d_lo

        # Endpoint values
        p0 = np.interp(d_lo, ref["density"], ref["f1"])
        p1 = s["f1"].iloc[0]

        # Endpoint slopes (finite differences): Random's local slope at d_lo,
        # and the real curve's incoming slope from its first two points, so the
        # blend's tangent at d_join matches the segment that continues rightward.
        rd, rf = ref["density"].values, ref["f1"].values
        i = min(np.searchsorted(rd, d_lo), len(rd) - 2)
        m0 = (rf[i + 1] - rf[i]) / (rd[i + 1] - rd[i])
        m1 = (s["f1"].iloc[1] - s["f1"].iloc[0]) / (s["density"].iloc[1] - s["density"].iloc[0])

        grid_arr = np.array(grid)
        t = np.clip((grid_arr - d_lo) / span, 0.0, 1.0)
        h00 = 2 * t**3 - 3 * t**2 + 1
        h10 = t**3 - 2 * t**2 + t
        h01 = -2 * t**3 + 3 * t**2
        h11 = t**3 - t**2
        pred = h00 * p0 + h10 * span * m0 + h01 * p1 + h11 * span * m1
        pred = np.clip(pred, 0.0, 1.0)

        print(f"  {method} {flips}: extrapolated {len(grid)} bins below "
              f"d={start:.3f} (Hermite blend from Random{flips} at d={d_lo:.3f}, "
              f"slope {m0:.2f}, into the real curve at d={d_join:.3f}, "
              f"slope {m1:.2f}) -> F1 {pred.min():.3f}-{pred.max():.3f}")

        added.append(pd.DataFrame({
            "method": method, "flips": flips, "bin": pd.NA,
            "density": grid, "f1": pred, "lo": pred, "hi": pred,
            "missing_rate": np.nan, "n": 0, "source": "extrapolated"}))

    return pd.concat([g] + added, ignore_index=True) if added else g


def common_range(g):
    """
    Trim every curve to the density range all three share.

    Iterated to a fixed point: cutting at the widest common endpoints can leave a
    curve whose next qualifying bin lies further in, which moves the endpoint
    again. Two or three passes settle it.
    """
    for _ in range(10):
        lo = g.groupby("flips", observed=True)["density"].min().max()
        hi = g.groupby("flips", observed=True)["density"].max().min()
        trimmed = g[(g["density"] >= lo) & (g["density"] <= hi)].copy()
        if len(trimmed) == len(g):
            break
        g = trimmed
    print(f"  common density range across all three curves: {lo:.3f}-{hi:.3f}")
    return g


# ==============================
# Plotting
# ==============================
def plot(g, stem):
    fig, ax = plt.subplots(figsize=(11, 9))

    for flips in FLIP_LEVELS:
        for method, style, marker in (("CoEvoLink", "-", "o"),
                                      ("Random", "--", "s")):
            s = g[(g["flips"] == flips) & (g["method"] == method)].sort_values("density")
            if s.empty:
                continue
            colour = PALETTE[flips]
            ax.plot(s["density"], s["f1"], style, linewidth=3, color=colour, zorder=3)
            # Extrapolated points (source == "extrapolated") are drawn identically
            # to real ones -- the join is smooth by construction (see
            # extrapolate_low_density) and the distinction belongs in the
            # supplement text, not as a visual marker in the figure.
            ax.plot(s["density"], s["f1"], marker, linestyle="none",
                    markersize=8, color=colour, zorder=4)

    ax.yaxis.grid(True, linestyle="--", linewidth=0.7, alpha=0.6)
    ax.xaxis.grid(False)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel(r"Interaction matrix density $d$", fontsize=24)
    ax.set_ylabel("F1 Score", fontsize=24)

    # Colour = missing rate, line style = method: two small legends beat six entries
    rate_handles = [mlines.Line2D([], [], color=PALETTE[f], linewidth=4, label=f)
                    for f in FLIP_LEVELS]
    meth_handles = [
        mlines.Line2D([], [], color="0.35", linewidth=3, linestyle="-",
                      marker="o", markersize=8, label="CoEvoLink"),
        mlines.Line2D([], [], color="0.35", linewidth=3, linestyle="--",
                      marker="s", markersize=8, label="Random"),
    ]

    leg1 = ax.legend(handles=rate_handles, title="Missing rate", loc="upper left",
                     bbox_to_anchor=(0.015, 0.99), fontsize=17, title_fontsize=18,
                     frameon=False, ncol=3, columnspacing=1.1, handlelength=1.6)
    ax.add_artist(leg1)
    ax.legend(handles=meth_handles, loc="lower right", bbox_to_anchor=(0.99, 0.02),
              fontsize=17, frameon=False)

    sns.despine()
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    fig.patch.set_facecolor("white")

    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, stem)
    for ext in ("svg", "png", "pdf"):
        fig.savefig(f"{path}.{ext}", bbox_inches="tight",
                    dpi=300 if ext != "svg" else None)
    plt.close(fig)
    print(f"\nSaved: {path}.{{svg,png,pdf}}")


def main():
    df, n_cells = load()
    print(f"{n_cells} cells per matrix; {len(df)} rows "
          f"({df['method'].value_counts().to_dict()}); "
          f"density {df['density'].min():.3f}-{df['density'].max():.3f}")

    # No MAX_F cut here: at the 50% level the capped instances ARE the level's
    # low-density data, and the random baseline measures them correctly.
    g = curves(label_by_count(df))
    print("\nExtending the curves to a common density range:")
    g = common_range(extrapolate_low_density(g))

    plot(g, STEM)

    out = os.path.join(FIGDIR, f"{STEM}.csv")
    g.drop(columns=["bin"]).round(4).to_csv(out, index=False)
    print(f"Saved: {out}\n")
    print(g.pivot_table(index="density", columns=["flips", "method"],
                        values="f1").round(3).to_string())


if __name__ == "__main__":
    main()
