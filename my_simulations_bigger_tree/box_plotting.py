import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ==============================
# Theme setup (Professor’s style)
# ==============================
plt.rcParams['xtick.labelsize'] = 22
plt.rcParams['ytick.labelsize'] = 22
plt.rcParams['axes.labelsize'] = 24

sns.set_context("notebook", font_scale=1.6)
sns.set_style("ticks")

# ==============================
# Load and merge
# ==============================
df_norm = pd.read_csv("simulation_all_metrics.csv")
df_elbow = pd.read_csv("simulation_all_metrics_elbow.csv")

df_elbow = df_elbow[~((df_elbow["precision"] < 0.5) | (df_elbow["f1"]  < 0.5))].copy()

df_norm["method"] = "Normal"
df_elbow["method"] = "Elbow"
df_all = pd.concat([df_norm, df_elbow], ignore_index=True)

# Only consider r01 = r10
df_eq = df_all[df_all["r01"] == df_all["r10"]].copy()
df_eq["rate"] = df_eq["r01"]
df_eq = df_eq[df_eq["rate"] <= 0.5]

# Keep only needed corruption levels
flip_map = {1000: "10%", 2000: "20%", 10000: "50%"}
df_eq = df_eq[df_eq["corrupt"].isin(flip_map.keys())]
df_eq["flips"] = df_eq["corrupt"].map(flip_map)

# ==============================
# Prepare output folder
# ==============================
os.makedirs("figures_boxplots_final", exist_ok=True)

# ==============================
# Colors for each metric (3 shades per metric)
# ==============================
palette_dict = {
    "precision": ["#aec7e8", "#1f77b4", "#08306b"],  # light→dark blue
    "recall": ["#aec7e8", "#1f77b4", "#08306b"],     # light→dark green
    "f1": ["#aec7e8", "#1f77b4", "#08306b"],         # light→dark red
}

# ==============================
# Plot each metric separately for Normal and Elbow
# ==============================
flip_levels = ["10%", "20%", "50%"]

for method in ["Normal", "Elbow"]:
    df_method = df_eq[df_eq["method"] == method]

    for metric in ["precision", "recall", "f1"]:
        fig, ax = plt.subplots(figsize=(9, 6))

        # Use hue for side-by-side grouping
        sns.boxplot(
            x="rate", y=metric, hue="flips",
            data=df_method,
            palette=palette_dict[metric],
            showfliers=False,
            ax=ax,
            width=0.6
        )

        # Optional: add jittered points for clarity
        sns.stripplot(
            x="rate", y=metric, hue="flips",
            data=df_method,
            dodge=True, alpha=0.35, size=5,
            palette=palette_dict[metric],
            ax=ax
        )

        # Clean legend and axes
        ax.yaxis.grid(True, linestyle="--", linewidth=0.7, alpha=0.6)
        ax.xaxis.grid(False)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel(r"Transition rate $\rho$", fontsize=24)
        ax.set_ylabel("Score", fontsize=24)

        # Remove duplicate legend (stripplot adds one more)
        handles, labels = ax.get_legend_handles_labels()
        n = len(flip_levels)
        ax.legend(
            handles[:n], flip_levels,
            title="Flip level",
            loc="upper center", bbox_to_anchor=(0.5, -0.12),
            ncol=3, fontsize=18, frameon=False, title_fontsize=20
        )

        sns.despine()
        plt.tight_layout(rect=[0, 0.05, 1, 1])
        fig.patch.set_facecolor("white")

        # Save figures
        fname_svg = f"figures_boxplots_final/{method.lower()}_{metric}.svg"
        fname_png = f"figures_boxplots_final/{method.lower()}_{metric}.png"
        plt.savefig(fname_svg, bbox_inches="tight")
        plt.savefig(fname_png, bbox_inches="tight", dpi=300)
        plt.close(fig)

        print(f"✅ Saved: {fname_svg} and {fname_png}")
