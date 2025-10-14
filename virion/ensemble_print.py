# produce_ensemble_curve.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json

import numpy as np
import pandas as pd
import json

def find_threshold_from_true_false(csv_path="WeightedEnsemble2020.csv"):
    """
    Finds approximate threshold separating TRUE and FALSE predictions
    based on the Ensemble.2 scores.
    """
    df = pd.read_csv(csv_path)

    # Clean and convert
    df["score"] = pd.to_numeric(df["Ensemble.2"], errors="coerce")
    df["label"] = (
        df["EnsembleBinary"]
        .astype(str)
        .str.strip()
        .str.lower()
        .replace({"true": True, "false": False})
    )

    # Extract TRUE and FALSE scores
    true_scores = df.loc[df["label"], "score"].dropna()
    false_scores = df.loc[~df["label"], "score"].dropna()

    if true_scores.empty or false_scores.empty:
        raise ValueError("No TRUE or FALSE labels found.")

    max_true = true_scores.max()
    min_false = false_scores.min()

    threshold = (max_true + min_false) / 2

    print(f"🔹 Highest TRUE score:  {max_true:.6f}")
    print(f"🔹 Lowest FALSE score:  {min_false:.6f}")
    print(f"✅ Approximate threshold ≈ {threshold:.6f}")


    return threshold

def check(
    ensemble_csv="WeightedEnsemble2020.csv",
    virion_json="virion_sets_named.json",
    model_name="ensemble"
):
    # --- Load data ---
    df = pd.read_csv(ensemble_csv, dtype=str)
    df["Ensemble.2"] = pd.to_numeric(df["Ensemble.2"], errors="coerce")
    with open(virion_json, "r", encoding="utf-8") as f:
        virion_sets = json.load(f)

    if model_name not in virion_sets:
        raise ValueError(f"Model '{model_name}' not found in {virion_json}")

    # --- Clean host names ---
    df["Sp_clean"] = df["Sp"].astype(str).str.strip().str.lower()

    # --- Virion host sets ---
    true_hosts = set(map(str.lower, virion_sets[model_name].get("true", [])))
    negative_hosts = set(map(str.lower, virion_sets[model_name].get("negative", [])))
    suspected_hosts = set(map(str.lower, virion_sets[model_name].get("suspected", [])))
    unlikely_hosts = set(map(str.lower, virion_sets[model_name].get("unlikely", [])))

    total_hosts = true_hosts | negative_hosts | suspected_hosts | unlikely_hosts

    df = df[df["Sp_clean"].isin(total_hosts)].copy()

        # Max score among true + suspected
    ts_scores = df[df["Sp_clean"].isin(true_hosts | suspected_hosts)]["Ensemble.2"]
    max_true_suspected = ts_scores.max() if not ts_scores.empty else None

    # Min score among negative + unlikely
    nu_scores = df[df["Sp_clean"].isin(negative_hosts | unlikely_hosts)]["Ensemble.2"]
    min_negative_unlikely = nu_scores.min() if not nu_scores.empty else None

    print(f"Max score (true + suspected): {max_true_suspected}")
    print(f"Min score (negative + unlikely): {min_negative_unlikely}")


def find_best_threshold_for_virion_sets(
    df, virion_json="virion_sets_named.json", score_col="Ensemble.2", model_name="ensemble"
):
    """
    Find the threshold (λ) where ensemble predicted positives best separate:
    - true, negative, suspected → predicted TRUE  (score < λ)
    - unlikely → predicted FALSE (score ≥ λ)
    
    Maximizes number of correctly classified hosts (not Jaccard).
    Uses binary search refinement around best coarse threshold.
    """

    # --- Load virion sets ---
    with open(virion_json, "r", encoding="utf-8") as f:
        virion_sets = json.load(f)

    if model_name not in virion_sets:
        raise ValueError(f"Model '{model_name}' not found in {virion_json}")

    true_hosts = set(map(str.lower, virion_sets[model_name].get("true", [])))
    negative_hosts = set(map(str.lower, virion_sets[model_name].get("negative", [])))
    suspected_hosts = set(map(str.lower, virion_sets[model_name].get("suspected", [])))
    unlikely_hosts = set(map(str.lower, virion_sets[model_name].get("unlikely", [])))

    positives = true_hosts | suspected_hosts
    negatives= negative_hosts | unlikely_hosts
    total_hosts = positives | negatives
    # --- Prepare ensemble data ---
    df = df.copy()
    if "Sp" not in df.columns:
        raise ValueError("DataFrame must contain column 'Sp' for host names.")
    if score_col not in df.columns:
        raise ValueError(f"DataFrame must contain column '{score_col}' for ensemble scores.")

    df["Sp_clean"] = df["Sp"].astype(str).str.strip().str.lower()
    df["score"] = df[score_col].astype(float)
    # Keep only hosts that are in virion total set
    df = df[df["Sp_clean"].isin(total_hosts)].dropna(subset=["score"])
    print(f"Loaded {len(df)} hosts overlapping with virion set.")

    # count how many are labeled True in the EnsembleBinary column
    true_count = (df["EnsembleBinary"] == "true").sum()
    print (true_count)

    # save df to csv
    df.to_csv("ensemble_virion_overlap.csv", index=False)

    # --- Define scoring function ---
    def correct_count(th):
        pred_true = set(df.loc[df["score"] < th, "Sp_clean"]) & positives
        pred_false = set(df.loc[df["score"] > th, "Sp_clean"]) & negatives
        return len(pred_true)

    # --- Coarse grid search ---
    # thresholds = np.arange(0, 1.001, 0.001)
    # best_thresh, best_correct = 0, -1
    # for t in thresholds:
    #     c = correct_count(t)
    #     if c > best_correct:
    #         best_thresh, best_correct = t, c

    # --- Binary search refinement ---
    low = 0
    high = 1
    eps = 1e-5

    while 1:
        mid = (low + high) / 2
        if correct_count(mid) > 134:
            high = mid
        elif correct_count(mid) < 134:
            low = mid
        else:
            low = mid
            break
    refined_thresh = low
    refined_correct = correct_count(refined_thresh)
    

    # --- Reporting ---
    print(f"⚙️  Refined correct: {refined_correct}")
    # # --- Save predictions ---
    # pred_true_final = set(df.loc[df["score"] <= refined_thresh, "Sp_clean"])
    # out_df = pd.DataFrame(
    #     sorted(pred_true_final),
    #     columns=["Predicted_Positive_Hosts"]
    # )
    # out_path = "ensemble_best_predpos_hosts.csv"
    # out_df.to_csv(out_path, index=False)
    # print(f"✅ Saved predicted positive hosts → {out_path}")

    # return refined_thresh



   
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
    # find_threshold_from_true_false(CSV_PATH)
    # find_best_threshold_for_virion_sets(df, virion_json="virion_sets_named.json")
    check(
    ensemble_csv="WeightedEnsemble2020.csv",
    virion_json="virion_sets_named.json",
    model_name="ensemble"
)