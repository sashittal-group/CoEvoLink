import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns   

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

# Filter: drop rows where both precision and recall == 0
# df = df[~((df["precision"] < 0.1 ) | (df["f1"] < 0.1) | (df["recall"] < 0.2))]

df.to_csv("simulation_all_metrics.csv", index=False)

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
# df = df[~((df["precision"] < 0.2 ) | (df["f1"] < 0.2) | (df["recall"] < 0.2))]
df = df[~((df["precision"] < 0.3 ))]

df.to_csv("simulation_all_metrics_elbow.csv", index=False)
print("Saved all_metrics.csv with shape:", df.shape)




# ==============================
# Theme setup (Professor’s style)
# ==============================
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 20

sns.set_context("notebook", font_scale=1.5)
sns.set_style("ticks")

# ==============================
# Load and merge
# ==============================
df_norm = pd.read_csv("simulation_all_metrics.csv")
df_elbow = pd.read_csv("simulation_all_metrics_elbow.csv")

df_norm["method"] = "Normal"
df_elbow["method"] = "Elbow"
df_all = pd.concat([df_norm, df_elbow], ignore_index=True)

# Only consider r01 = r10
df_eq = df_all[df_all["r01"] == df_all["r10"]].copy()
df_eq["rate"] = df_eq["r01"]

# ==============================
# Data reshape for boxplot
# ==============================
df_melted = df_eq.melt(
    id_vars=["method", "corrupt", "rate"],
    value_vars=["precision", "recall", "f1"],
    var_name="metric",
    value_name="score"
)

# ==============================
# Output folder
# ==============================
os.makedirs("figures_boxplots", exist_ok=True)

# ==============================
# Plot loop
# ==============================
palette_metrics = {
    "precision": "#1f77b4",  # blue
    "recall": "#2ca02c",     # green
    "f1": "#d62728",         # red
}

for method in ["Normal", "Elbow"]:
    df_method = df_melted[df_melted["method"] == method]
    flip_levels = sorted(df_method["corrupt"].unique())

    for corrupt in flip_levels:
        df_plot = df_method[df_method["corrupt"] == corrupt]

        fig, ax = plt.subplots(figsize=(9, 6))

        sns.boxplot(
            x="rate", y="score", hue="metric",
            data=df_plot, palette=palette_metrics,
            showfliers=False, ax=ax
        )

        sns.stripplot(
            x="rate", y="score", hue="metric",
            data=df_plot, dodge=True, alpha=0.4,
            linewidth=1, jitter=0.1,
            palette=palette_metrics, ax=ax
        )

        # Light horizontal grid only
        ax.yaxis.grid(True, linestyle="--", linewidth=0.7, alpha=0.6)
        ax.xaxis.grid(False)

        # Clean legend and move it to the bottom
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles[0:3], ['Precision', 'Recall', 'F1'],
            loc='upper center', bbox_to_anchor=(0.5, -0.15),
            ncol=3, fontsize=14, frameon=False
        )

        # Axis and title settings
        ax.set_ylim(0, 1.05)
        ax.set_xlabel(r"r$_{01}$ = r$_{10}$", fontsize=20)
        ax.set_ylabel("Score", fontsize=20)

        if corrupt == 500:
            percentage = 2
        elif corrupt == 1000:
            percentage = 5
        elif corrupt == 5000:
            percentage = 20
        elif corrupt == 10000:
            percentage = 50
        else:
            percentage = corrupt

        ax.set_title(f"{percentage}% flips", fontsize=22, fontweight="bold")

        sns.despine()
        plt.tight_layout(rect=[0, 0.05, 1, 1])  # leave space at bottom for legend
        fig.patch.set_facecolor("white")

        # Save plots
        fname_svg = f"figures_boxplots/{method.lower()}_corrupt{corrupt}.svg"
        fname_png = f"figures_boxplots/{method.lower()}_corrupt{corrupt}.png"
        plt.savefig(fname_svg, bbox_inches="tight")
        plt.savefig(fname_png, bbox_inches="tight", dpi=300)
        plt.close(fig)

        print(f"✅ Saved: {fname_svg} and {fname_png}")

