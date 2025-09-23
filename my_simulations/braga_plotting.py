#!/usr/bin/env python3
# plotting_braga.py
import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = "results/braga"
OUT_CSV = "braga_all_metrics.csv"

records = []

# Walk through results/braga/
for root, dirs, files in os.walk(BASE_DIR):
    if "metrics.json" not in files:
        continue

    metrics_path = os.path.join(root, "metrics.json")

    # path pattern: results/braga/{seed}/b{beta}_c{clock}/metrics.json
    parts = os.path.normpath(root).split(os.sep)
    if len(parts) < 4:
        continue
    try:
        seed = int(parts[2])  # {seed}
    except Exception:
        continue

    subdir = parts[3]  # "b{beta}_c{clock}"
    if not subdir.startswith("b") or "_c" not in subdir:
        continue

    try:
        beta_str, clock_str = subdir.split("_c")
        beta = int(beta_str.replace("b", ""))
        clock = float(clock_str)
    except Exception as e:
        print(f"[WARN] Could not parse beta/clock from {subdir}: {e}")
        continue

    try:
        with open(metrics_path) as fh:
            metrics = json.load(fh)
    except Exception as e:
        print(f"[WARN] Failed to read {metrics_path}: {e}")
        continue

    if not all(k in metrics for k in ("precision", "recall", "f1")):
        print(f"[WARN] Incomplete metrics in {metrics_path}; skipping")
        continue

    record = {
        "seed": seed,
        "beta": beta,
        "clock": clock,
        "precision": float(metrics.get("precision", 0.0)),
        "recall": float(metrics.get("recall", 0.0)),
        "f1": float(metrics.get("f1", 0.0)),
        "TP": int(metrics.get("TP", 0)),
        "FP": int(metrics.get("FP", 0)),
        "FN": int(metrics.get("FN", 0)),
    }
    records.append(record)

# Build dataframe
df = pd.DataFrame(records)
if df.empty:
    raise SystemExit("No BRAGA metrics found in results/braga/ — nothing to plot.")

# Drop rows with both precision and recall == 0 (failed runs)
df = df[~((df["precision"] < 0.35) & (df["recall"]  < 0.35))].copy()
df.to_csv(OUT_CSV, index=False)
print(f"Saved {OUT_CSV} with shape: {df.shape}")

# Average across seeds
df_avg = df.groupby(["beta", "clock"])[["precision", "recall", "f1"]].mean().reset_index()

# Plot heatmaps for each metric
for metric in ["precision", "recall", "f1"]:
    pivot = df_avg.pivot(index="beta", columns="clock", values=metric)
    plt.figure(figsize=(8, 6))
    sns.heatmap(
        pivot,
        annot=True, fmt=".2f", cmap="viridis", vmin=0, vmax=1,
        cbar_kws={"label": metric.capitalize()}
    )
    plt.title(f"BRAGA {metric.capitalize()} (avg over seeds)")
    plt.ylabel("Beta")
    plt.xlabel("Clock")
    plt.tight_layout()
    fname = f"braga_{metric}.png"
    plt.savefig(fname)
    plt.close()
    print(f"Wrote {fname}")

print("Done.")
