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
df = df[~((df["precision"] < 0.3) | (df["recall"] < 0.3) | (df["f1"] < 0.3))]

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
df = df[~((df["precision"] < 0.3 ) & (df["f1"] < 0.3))]

df.to_csv("simulation_all_metrics_elbow.csv", index=False)
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
    plt.savefig(f"{metric}_by_corrupt_elbow.svg")
    plt.close()