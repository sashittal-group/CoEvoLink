# produce_ensemble_curve.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Path to the CSV (change if needed)
CSV_PATH = "WeightedEnsemble2020.csv"

def compute_predpos_acc(df, score_col="Ensemble.2", truth_col="Binary", thresholds=None):
    """
    Compute PredPos and Acc for the supplied thresholds.
    - df: DataFrame
    - score_col: name of numeric score column (Ensemble.2)
    - truth_col: name of binary truth column (1/0), used to compute Acc
    - thresholds: 1D array of thresholds; predicted positive = score < thresh
    Returns: DataFrame with columns thresh, predpos, acc, count_pred, count_acc
    """
    if thresholds is None:
        thresholds = np.arange(0, 1001) / 1000.0  # default: ensemble thresholds (0.000 .. 1.000), length 1001

    scores = pd.to_numeric(df[score_col], errors="coerce")
    valid_mask = scores.notna()
    scores = scores[valid_mask].values
    truths = pd.to_numeric(df.loc[valid_mask, truth_col], errors="coerce").fillna(0).astype(int).values

    rows = []
    counts = []
    acc_counts = []
    for t in thresholds:
        pred_mask = (scores < t)             # R used "< thresh"
        count_pred = int(pred_mask.sum())
        count_acc = int((pred_mask & (truths == 1)).sum())
        counts.append(count_pred)
        acc_counts.append(count_acc)
        rows.append({"thresh": float(t),
                     "count_pred": count_pred,
                     "count_acc": count_acc})
    counts = np.array(counts, dtype=float)
    acc_counts = np.array(acc_counts, dtype=float)

    # Normalize exactly like R: divide by max(counts) and max(acc_counts)
    max_count = counts.max() if counts.max() > 0 else 1.0
    max_acc = acc_counts.max() if acc_counts.max() > 0 else 1.0
    predpos = counts / max_count
    acc = acc_counts / max_acc

    out_rows = []
    for i, t in enumerate(thresholds):
        out_rows.append({
            "thresh": float(t),
            "predpos": float(predpos[i]),
            "acc": float(acc[i]),
            "count_pred": int(counts[i]),
            "count_acc": int(acc_counts[i])
        })
    return pd.DataFrame(out_rows)

def plot_curve(df_curve, savepath="ensemble_predpos_vs_acc.png"):
    plt.figure(figsize=(7,6))
    plt.plot(df_curve["predpos"] * 100, df_curve["acc"] * 100, lw=1.5, color="black")
    plt.xlabel("Predicted positivity rate (%)")
    plt.ylabel("New hosts correctly identified (%)")
    plt.title("Ensemble: PredPos vs Acc (matches R/X_BetterEnsemble.R)")
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(savepath, dpi=300)
    print("Saved plot to", savepath)
    plt.close()

if __name__ == "__main__":
    df = pd.read_csv(CSV_PATH, dtype=str).fillna("")
    if "Ensemble.2" not in df.columns:
        raise SystemExit("CSV must contain column 'Ensemble.2' with numeric ensemble scores.")
    if "Binary" not in df.columns:
        print("Warning: 'Binary' column not present. Acc will be zero if no positives found.")

    # EXACT thresholds used by the R ensemble run:
    thresholds_ensemble = np.arange(0, 10001) / 10000.0  # 0.000, 0.001, ..., 1.000 (length 1001)

    # If you later need to compute per-model thresholds (R used an extra point), use:
    # thresholds_models = np.arange(0, 1002) / 1000.0  # 0.000 .. 1.001 (length 1002)

    curve_df = compute_predpos_acc(df, score_col="Ensemble.2", truth_col="Binary", thresholds=thresholds_ensemble)
    curve_df.to_csv("ensemble_predpos_acc_curve.csv", index=False)
    plot_curve(curve_df, savepath="ensemble_predpos_vs_acc.png")

    # Print sample and find the threshold maximizing (Acc - PredPos) like R's tdf$Diff
    curve_df["diff"] = curve_df["acc"] - curve_df["predpos"]
    best = curve_df.loc[curve_df["diff"].idxmax()]
    print("\nBest threshold (max Acc - PredPos):")
    print(best)
    print("\nSaved curve CSV -> ensemble_predpos_acc_curve.csv")