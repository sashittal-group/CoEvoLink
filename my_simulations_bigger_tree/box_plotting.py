import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns   

base_dir = "results"
records = []

# ---- Collect metrics.json ----
for root, dirs, files in os.walk(base_dir):
    if "metrics.json" in files:
        metrics_path = os.path.join(root, "metrics.json")
        with open(metrics_path) as f:
            metrics = json.load(f)
        
        parts = root.split(os.sep)
        if parts[1] == 'braga':  # Skip braga results
            continue

        seed = int(parts[1])
        corrupt = int(parts[2].replace("corrupt", ""))
        rate_part = parts[3]  # "r01X_r10Y"
        r01 = float(rate_part.split("_")[0].replace("r01", ""))
        r10 = float(rate_part.split("_")[1].replace("r10", ""))

        record = {
            "seed": seed,
            "corrupt": corrupt,
            "r01": r01,
            "r10": r10,
            **metrics
        }
        records.append(record)

df = pd.DataFrame(records)

# ---- Filter out poor runs ----
# df = df[~((df["precision"] < 0.0) & (df["recall"] < 0.0))]
df.to_csv("simulation_all_metrics.csv", index=False)
print("Saved all_metrics.csv with shape:", df.shape)

# ---- Add rate pair label ----
df["rate_pair"] = df.apply(lambda row: f"({row.r01},{row.r10})", axis=1)

# ---- Filter only r01 > r10 ----
df_asym = df[df["r01"] > df["r10"]].copy()

# ---- Melt for plotting ----
df_long = df_asym.melt(
    id_vars=["corrupt", "rate_pair", "seed"],
    value_vars=["precision", "recall", "f1"],
    var_name="metric",
    value_name="score"
)

# ---- Filter for specific corruption (e.g. 1000) ----
df_long = df_long[df_long["corrupt"] == 1000]

# ---- Single boxplot for corrupt = 1000 ----
g = sns.catplot(
    data=df_long,
    x="rate_pair", y="score",
    hue="metric",
    kind="box",
    dodge=0.8,
    width=0.7,
    showfliers=False,
    height=7, aspect=1
)

# ---- Formatting ----
for ax in g.axes.flat:
    ax.tick_params(axis="x", rotation=45)

g.set_axis_labels("rate (r01, r10)", "Score")
g.set_titles("Random Flipping = 1000")
g.set(ylim=(0, 1))
g.fig.subplots_adjust(top=0.85)
g.fig.suptitle("Precision / Recall / F1 across (r01, r10) pairs", fontsize=14)

plt.savefig("metrics_boxplot.png", dpi=300, bbox_inches="tight")
plt.close()

print("✅ Saved metrics_boxplot.png")

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

# Filter: drop rows where both precision and recall < 0.5
# df = df[~((df["precision"] < 0.0) & (df["recall"] < 0.0))]

df.to_csv("simulation_all_metrics_elbow.csv", index=False)
print("Saved all_metrics.csv with shape:", df.shape)

# ---- Add rate pair label ----
df["rate_pair"] = df.apply(lambda row: f"({row.r01},{row.r10})", axis=1)

# ---- Filter only r01 > r10 ----
df_asym = df[df["r01"] > df["r10"]].copy()

# ---- Melt for plotting ----
df_long = df_asym.melt(
    id_vars=["corrupt", "rate_pair", "seed"],
    value_vars=["precision", "recall", "f1"],
    var_name="metric",
    value_name="score"
)

# ---- Filter for specific corruption (e.g. 1000) ----
df_long = df_long[df_long["corrupt"] == 1000]

# ---- Single boxplot for corrupt = 1000 ----
g = sns.catplot(
    data=df_long,
    x="rate_pair", y="score",
    hue="metric",
    kind="box",
    dodge=0.8,
    width=0.7,
    showfliers=False,
    height=7, aspect=1
)

# ---- Formatting ----
for ax in g.axes.flat:
    ax.tick_params(axis="x", rotation=45)

g.set_axis_labels("rate (r01, r10)", "Score")
g.set_titles("Random Flipping = 1000")
g.set(ylim=(0, 1))
g.fig.subplots_adjust(top=0.85)
g.fig.suptitle("Precision / Recall / F1 across (r01, r10) pairs", fontsize=14)

plt.savefig("metrics_boxplot_elbow.svg", dpi=300, bbox_inches="tight")
plt.close()

print("✅ Saved metrics_boxplot_elbow.svg")