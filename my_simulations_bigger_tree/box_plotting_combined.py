import os
import json
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


base_dir = "results"
records = []

# Walk through results recursively
for root, dirs, files in os.walk(base_dir):
    if "metrics.json" in files:
        metrics_path = os.path.join(root, "metrics.json")
        with open(metrics_path) as f:
            metrics = json.load(f)
        
        # Parse path: results/{seed}/corrupt{corrupt}/r01{r01}_r10{r10}/
        parts = root.split(os.sep)
        if (parts[1] == 'braga'):  # Skip braga results'
            continue
        seed = int(parts[1])  # "results/{seed}"
        corrupt = int(parts[2].replace("corrupt", ""))
        rate_part = parts[3]  # "r01X_r10Y"
        r01 = float(rate_part.split("_")[0].replace("r01", ""))
        r10 = float(rate_part.split("_")[1].replace("r10", ""))

        record = {
            "seed": seed,
            "corrupt": corrupt,
            "r01": r01,
            "r10": r10,
            **metrics  # precision, recall, f1, TP, FP, FN
        }
        records.append(record)

df = pd.DataFrame(records)


df.to_csv("simulation_all_metrics.csv", index=False)
print("Saved all_metrics.csv with shape:", df.shape)

# Add rate pair label
df["rate_pair"] = df.apply(lambda row: f"({row.r01},{row.r10})", axis=1)

# Average only over filtered data
df_avg = df.groupby(["corrupt", "rate_pair"])[["precision", "recall", "f1"]].mean().reset_index()

# Plot
for metric in ["precision", "recall", "f1"]:
    plt.figure(figsize=(10,6))
    sns.barplot(data=df_avg, x="rate_pair", y=metric, hue="corrupt", ci=None)
    plt.title(f"{metric.capitalize()} across corruption levels")
    plt.ylim(0,1)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{metric}_by_corrupt.svg")
    plt.close()


base_dir = "results"
records = []

# Walk through results recursively
for root, dirs, files in os.walk(base_dir):
    if "elbow_metrics.json" in files:
        metrics_path = os.path.join(root, "elbow_metrics.json")
        with open(metrics_path) as f:
            metrics = json.load(f)
        
        # Parse path: results/{seed}/corrupt{corrupt}/r01{r01}_r10{r10}/
        parts = root.split(os.sep)
        if (parts[1] == 'braga'):  # Skip braga results'
            continue
        seed = int(parts[1])  # "results/{seed}"
        corrupt = int(parts[2].replace("corrupt", ""))
        rate_part = parts[3]  # "r01X_r10Y"
        r01 = float(rate_part.split("_")[0].replace("r01", ""))
        r10 = float(rate_part.split("_")[1].replace("r10", ""))

        record = {
            "seed": seed,
            "corrupt": corrupt,
            "r01": r01,
            "r10": r10,
            **metrics  # precision, recall, f1, TP, FP, FN
        }
        records.append(record)

df = pd.DataFrame(records)

# Filter: drop rows where both precision and recall == 0
df = df[~((df["precision"] < 0.35 ) & (df["f1"] < 0.4))]

df.to_csv("simulation_all_metrics_elbow.csv", index=False)
print("Saved all_metrics_elbow.csv with shape:", df.shape)

# ==============================
# Theme setup
# ==============================
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 22

sns.set_context("notebook", font_scale=1.5)
sns.set_style("ticks")

# ==============================
# Load and merge
# ==============================
df_norm = pd.read_csv("simulation_all_metrics.csv")
df_elbow = pd.read_csv("simulation_all_metrics_elbow.csv")

# df_elbow = df_elbow[~((df_elbow["precision"] < 0.3) | (df_elbow["f1"]  < 0.3))].copy()

df_norm["method"] = "Normal"
df_elbow["method"] = "Elbow"
df_all = pd.concat([df_norm, df_elbow], ignore_index=True)

# Only consider r01 = r10
df_eq = df_all[df_all["r01"] == df_all["r10"]].copy()
df_eq["rate"] = df_eq["r01"]
df_eq = df_eq[df_eq["rate"] <= 0.5]

# Keep only desired corruption levels
flip_map = {1000: "10%", 2000: "20%", 10000: "50%"}
df_eq = df_eq[df_eq["corrupt"].isin(flip_map.keys())]
df_eq["flips"] = df_eq["corrupt"].map(flip_map)

# ==============================
# Output folder
# ==============================
os.makedirs("figures_boxplots_combined", exist_ok=True)

# ==============================
# Shared color palettes
# ==============================
palette_elbow = ["#aec7e8", "#1f77b4", "#08306b"]  # light→dark blue
palette_normal = ["#fcae91", "#fb6a4a", "#a50f15"]   # light→dark red

flip_levels = ["10%", "20%", "50%"]
metrics = ["precision", "recall", "f1"]
metric_labels = {"precision": "Precision", "recall": "Recall", "f1": "F1-score"}

# ==============================
# Combined figure per method
# ==============================
for method, palette in zip(["Normal", "Elbow"], [palette_normal, palette_elbow]):
    df_method = df_eq[df_eq["method"] == method]

    # Remove sharey=True
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    y_min, y_max = 0, 1.05  # same scale for all subplots

    for ax, metric in zip(axes, metrics):
        # --- Boxplot ---
        sns.boxplot(
            x="rate", y=metric, hue="flips",
            data=df_method,
            palette=palette,
            showfliers=False,
            ax=ax,
            width=0.6
        )

        # --- Points (stripplot) ---
        sns.stripplot(
            x="rate", y=metric, hue="flips",
            data=df_method,
            dodge=True,
            alpha=0.35, size=5,
            palette=palette,
            ax=ax
        )

        # --- Axis & Style ---
        ax.yaxis.grid(True, linestyle="--", linewidth=0.7, alpha=0.6)
        ax.xaxis.grid(False)
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel(r"Transition rate $\rho$", fontsize=22)
        ax.set_ylabel(metric_labels[metric], fontsize=22)
        ax.get_legend().remove()  # Remove duplicate legends

    # --- Shared legend at bottom ---
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles[:len(flip_levels)], flip_levels,
        title="Flip level",
        loc="lower center",
        ncol=3,
        fontsize=18,
        title_fontsize=20,
        frameon=False,
        bbox_to_anchor=(0.5, -0.15)
    )

    sns.despine()
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.15, top=0.95, wspace=0.28)
    fig.patch.set_facecolor("white")

    # --- Save figure ---
    fname_svg = f"figures_boxplots_combined/{method.lower()}_combined.svg"
    fname_png = f"figures_boxplots_combined/{method.lower()}_combined.png"
    plt.savefig(fname_svg, bbox_inches="tight")
    plt.savefig(fname_png, bbox_inches="tight", dpi=300)
    plt.close(fig)

    print(f"✅ Saved: {fname_svg} and {fname_png}")
