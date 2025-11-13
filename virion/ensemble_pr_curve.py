#!/usr/bin/env python3
"""
Compute precision-recall curve for ensemble predictions using virion_sets_named.json mapping.

Rules followed:
- ground positives = union of 'true' and 'negative' lists from the JSON under key 'ensemble'
- ground negatives = union of 'suspected' and 'unlikely' lists from the JSON under key 'ensemble'
- only species present in the JSON 'ensemble' lists are considered
- a predicted positive is a row where the ensemble score (default 'Ensemble.2') < threshold (matches the provided script)

Outputs:
- CSV: ensemble_pr_curve.csv (thresh, precision, recall, TP, FP, FN, n_predicted, n_ground_pos)
- PNG: ensemble_precision_recall.png
"""
import argparse
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_json_sets(json_path: str):
    with open(json_path, 'r') as f:
        data = json.load(f)
    if 'ensemble' not in data:
        raise KeyError("JSON must contain top-level key 'ensemble'")
    ens = data['ensemble']
    # helper to lower and strip
    def norm_list(l):
        return [s.strip().lower() for s in l]

    true = set(norm_list(ens.get('true', [])))
    negative = set(norm_list(ens.get('negative', [])))
    suspected = set(norm_list(ens.get('suspected', [])))
    unlikely = set(norm_list(ens.get('unlikely', [])))

    ground_pos = true.union(negative)
    ground_neg = suspected.union(unlikely)
    ensemble_hosts = true.union(negative).union(suspected).union(unlikely)
    return {
        'true': true,
        'negative': negative,
        'suspected': suspected,
        'unlikely': unlikely,
        'ground_pos': ground_pos,
        'ground_neg': ground_neg,
        'ensemble_hosts': ensemble_hosts,
    }


def compute_pr_curve(df: pd.DataFrame, score_col: str, host_col: str, sets: dict, thresholds: np.ndarray):
    # normalize species column
    df = df.copy()
    df['__host_norm'] = df[host_col].astype(str).str.strip().str.lower()

    # filter to ensemble hosts only
    in_ensemble = df['__host_norm'].isin(sets['ensemble_hosts'])
    df = df.loc[in_ensemble].reset_index(drop=True)

    scores = pd.to_numeric(df[score_col], errors='coerce').values
    hosts = df['__host_norm'].values

    rows = []
    # Precompute masks for membership
    is_ground_pos = np.array([h in sets['ground_pos'] for h in hosts])
    is_ground_neg = np.array([h in sets['ground_neg'] for h in hosts])

    for t in thresholds:
        pred_pos = scores < float(t)
        TP = int(((pred_pos) & (is_ground_pos)).sum())
        FP = int(((pred_pos) & (is_ground_neg)).sum())
        FN = int(((~pred_pos) & (is_ground_pos)).sum())
        TN = int(((~pred_pos) & (is_ground_neg)).sum())
        n_pred = int(pred_pos.sum())
        n_ground_pos = int(is_ground_pos.sum())

        precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
        recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0

        rows.append({
            'thresh': float(t),
            'precision': float(precision),
            'recall': float(recall),
            'TP': TP,
            'FP': FP,
            'FN': FN,
            'TN': TN,
            'n_predicted': n_pred,
            'n_ground_pos': n_ground_pos,
        })

    return pd.DataFrame(rows)


def plot_pr(df_pr: pd.DataFrame, outpath: str):
    plt.figure(figsize=(7,6))
    plt.plot(df_pr['recall'], df_pr['precision'], lw=1.5, color='blue')
    plt.scatter(df_pr['recall'], df_pr['precision'], s=8, alpha=0.5)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Ensemble: Precision-Recall curve')
    plt.grid(True, alpha=0.3)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    print('Saved PR plot to', outpath)
    plt.close()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--csv', default='WeightedEnsemble2020.csv', help='CSV with columns for host and ensemble score')
    p.add_argument('--json', default='virion_sets_named.json', help='JSON file with ensemble sets')
    p.add_argument('--host-col', default='Sp', help='Column name for species/host names in CSV')
    p.add_argument('--score-col', default='Ensemble.2', help='Column name for numeric ensemble score')
    p.add_argument('--out-csv', default='ensemble_pr_curve.csv')
    p.add_argument('--out-plot', default='ensemble_precision_recall.png')
    p.add_argument('--threshold-step', type=float, default=1/1000.0, help='Step for thresholds (default 0.0001)')
    args = p.parse_args()

    csv_path = Path(args.csv)
    json_path = Path(args.json)

    if not csv_path.exists():
        raise SystemExit(f'CSV not found: {csv_path}')
    if not json_path.exists():
        raise SystemExit(f'JSON not found: {json_path}')

    df = pd.read_csv(csv_path, dtype=str).fillna('')
    if args.host_col not in df.columns:
        raise SystemExit(f"Host column '{args.host_col}' not found in CSV columns: {list(df.columns)}")
    if args.score_col not in df.columns:
        raise SystemExit(f"Score column '{args.score_col}' not found in CSV columns: {list(df.columns)}")

    sets = load_json_sets(str(json_path))

    thresholds = np.arange(0, 1.0 + args.threshold_step/2, args.threshold_step)

    df_pr = compute_pr_curve(df, score_col=args.score_col, host_col=args.host_col, sets=sets, thresholds=thresholds)

    df_pr.to_csv(args.out_csv, index=False)
    print('Saved PR CSV to', args.out_csv)

    plot_pr(df_pr, args.out_plot)

    # print best F1
    f1 = (2 * df_pr['precision'] * df_pr['recall']) / (df_pr['precision'] + df_pr['recall']).replace({0: np.nan})
    best_idx = f1.idxmax()
    if pd.notna(best_idx):
        print('\nBest threshold by F1:')
        print(df_pr.loc[int(best_idx)].to_dict())


if __name__ == '__main__':
    main()
