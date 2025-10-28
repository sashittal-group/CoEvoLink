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
df = df[~((df["precision"] == 0) & (df["recall"] == 0))]

df.to_csv("simulation_all_metrics.csv", index=False)
print("Saved all_metrics.csv with shape:", df.shape)

# Add rate pair label
df["rate_pair"] = df.apply(lambda row: f"({row.r01},{row.r10})", axis=1)

# Average only over filtered data
df_avg = df.groupby(["corrupt", "rate_pair"])[["precision", "recall", "f1"]].mean().reset_index()

# Keep only symmetric r01=r10
df_sym = df_avg[df_avg["rate_pair"].apply(lambda x: x[1:-1].split(",")[0] == x[1:-1].split(",")[1])].copy()
# Keep only symmetric r01=r10
df_sym = df[df["r01"] == df["r10"]].copy()
df_sym["r"] = df_sym["r01"]

# Melt into long format for plotting (do NOT average, keep per-seed data!)
df_long = df_sym.melt(
    id_vars=["corrupt", "r", "seed"],
    value_vars=["precision", "recall", "f1"],
    var_name="metric",
    value_name="score"
)

# Faceted boxplots, one panel per corruption
g = sns.catplot(
    data=df_long,
    x="r", y="score",
    hue="metric",
    col="corrupt",
    kind="box",
    dodge=True,
    width=0.7,          # narrower boxes → spacing between them
    gap=0.4,            # add gap between boxes within each group
    showfliers=False,    # remove "o" outliers
    height=7, aspect=1
)

# Formatting
g.set_axis_labels("rate 01 = rate 10", "Score")
g.set_titles("Random Flipping = {col_name}")
g.set(ylim=(0,1))
g.fig.subplots_adjust(top=0.85, wspace=0.15)  # wspace → spacing between panels
g.fig.suptitle("Precision / Recall / F1 across rates and flipping", fontsize=14)

plt.savefig("metrics_boxplot.png", dpi=300, bbox_inches="tight")
plt.close()

# Filter for corrupt = 100
df_long = df_long[df_long["corrupt"] == 100]

# Single boxplot for corrupt = 100
g = sns.catplot(
    data=df_long,
    x="r", y="score",
    hue="metric",
    kind="box",
    dodge=0.8,
    width=0.7,          # narrower boxes → spacing between them
    showfliers=False,    # remove "o" outliers
    height=7, aspect=1
)

# Formatting
g.set_axis_labels("rate 01 = rate 10", "Score")
g.set_titles("Random Flipping = 100")
g.set(ylim=(0,1))
g.fig.subplots_adjust(top=0.85)
g.fig.suptitle("Precision / Recall / F1 across rates for flipping = 100", fontsize=14)

plt.savefig("metrics_boxplot.svg", dpi=300, bbox_inches="tight")
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
df = df[~((df["precision"] < 0.3) & (df["f1"] < 0.3))]

df.to_csv("simulation_all_metrics_elbow.csv", index=False)
print("Saved all_metrics.csv with shape:", df.shape)

# Add rate pair label
df["rate_pair"] = df.apply(lambda row: f"({row.r01},{row.r10})", axis=1)

# Average only over filtered data
df_avg = df.groupby(["corrupt", "rate_pair"])[["precision", "recall", "f1"]].mean().reset_index()


# Keep only symmetric r01=r10
df_sym = df_avg[df_avg["rate_pair"].apply(lambda x: x[1:-1].split(",")[0] == x[1:-1].split(",")[1])].copy()
# Keep only symmetric r01=r10
df_sym = df[df["r01"] == df["r10"]].copy()
df_sym["r"] = df_sym["r01"]

# Melt into long format for plotting (do NOT average, keep per-seed data!)
df_long = df_sym.melt(
    id_vars=["corrupt", "r", "seed"],
    value_vars=["precision", "recall", "f1"],
    var_name="metric",
    value_name="score"
)

# Filter for corrupt = 100
df_long = df_long[df_long["corrupt"] == 100]

# Single boxplot for corrupt = 100
g = sns.catplot(
    data=df_long,
    x="r", y="score",
    hue="metric",
    kind="box",
    dodge=0.8,
    width=0.7,          # narrower boxes → spacing between them
    showfliers=False,    # remove "o" outliers
    height=7, aspect=1
)

# Formatting
g.set_axis_labels("rate 01 = rate 10", "Score")
g.set_titles("Random Flipping = 100")
g.set(ylim=(0,1))
g.fig.subplots_adjust(top=0.85)
g.fig.suptitle("Precision / Recall / F1 across rates for flipping = 100", fontsize=14)

plt.savefig("metrics_boxplot_elbow.png", dpi=300, bbox_inches="tight")
plt.close()