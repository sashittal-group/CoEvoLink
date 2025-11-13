import json
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
#  Load data
# -----------------------------
methods = ["pblks", "phist", "wish"]
datasets = ["cow feces", "human gut", "wastewater"]

# your colors
colors = {
    "pblks": "#ff9999",   # light red (PB-LKS)
    "wish": "#d66ce7f6",  # light pink (WISH)
    "phist": "#5aca83",   # light green (PHIST)
    "coevo": "#1f77b4",   # dark blue (CoEvoLink)
}

data = {}
for method in methods:
    with open(f"{method}_top1.json") as f:
        data[f"{method}_top1"] = json.load(f)
    with open(f"{method}_top3.json") as f:
        data[f"{method}_top3"] = json.load(f)
    with open(f"coevo_top1_{method}.json") as f:
        data[f"coevo_top1_{method}"] = json.load(f)
    with open(f"coevo_top3_{method}.json") as f:
        data[f"coevo_top3_{method}"] = json.load(f)

# -----------------------------
#  Plot settings
# -----------------------------
bar_width = 0.18
x = np.arange(len(datasets)) * 1.2

plt.rcParams.update({"font.size": 22})

# ↓ shorter figure: reduce height from 8 → 5
fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)

offset_1 = -bar_width * 1.4
offset_3 = bar_width * 1.4

titles = {
    "pblks": "PB-LKS",
    "wish": "WISH",
    "phist": "PHIST"
}

# -----------------------------
#  Loop over methods
# -----------------------------
for ax, method in zip(axes, methods):
    top1_ours = data[f"{method}_top1"]
    top3_ours = data[f"{method}_top3"]
    top1_coevo = data["coevo_top1_" + method]
    top3_coevo = data["coevo_top3_" + method]

    # Bars
    bars1_ours = ax.bar(x + offset_1 - bar_width / 2, top1_ours, width=bar_width,
                        color=colors[method], label=titles[method])
    bars1_coevo = ax.bar(x + offset_1 + bar_width / 2, top1_coevo, width=bar_width,
                         color=colors["coevo"], label="CoEvoLink")

    bars3_ours = ax.bar(x + offset_3 - bar_width / 2, top3_ours, width=bar_width,
                        color=colors[method])
    bars3_coevo = ax.bar(x + offset_3 + bar_width / 2, top3_coevo, width=bar_width,
                         color=colors["coevo"])

    # -----------------------------
    # Add accuracy labels on top
    # -----------------------------
    def add_labels(bars):
        for bar in bars:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                height + 0.02,
                f"{height:.2f}",
                ha="center", va="bottom",
                rotation=90,           # vertical orientation
                fontsize=18,
                # fontweight="bold"      # optional: makes it clearer
            )


    for group in [bars1_ours, bars1_coevo, bars3_ours, bars3_coevo]:
        add_labels(group)

    # Labels and formatting
    ax.set_xticks(x)
    ax.set_ylim(0, 1.1)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.legend(fontsize=18, ncol=2, loc="upper center")

    # Dataset + top1/top3 labels
    for i, dataset in enumerate(datasets):
        x1 = x[i] + offset_1
        x3 = x[i] + offset_3
        y_offset = -0.08
        ax.text(x[i], y_offset-0.03, dataset.capitalize(), ha="center", va="top", fontsize=22)
        ax.text(x1, y_offset - 0.01, "top-1", ha="center", va="bottom", fontsize=21)
        ax.text(x3, y_offset - 0.01, "top-3", ha="center", va="bottom", fontsize=21)

    ax.set_xticklabels([])
    if ax is axes[0]:
        ax.set_ylabel("Accuracy", fontsize=22)
    else:
        ax.set_ylabel("")

# Layout
plt.subplots_adjust(wspace=0.1, top=0.85, bottom=0.25)
plt.tight_layout(rect=[0, 0, 1, 0.9])

plt.savefig("combined_methods_top1_top3_accuracy.pdf", dpi=300)
plt.close()
print("✅ Figure saved with compact height and clean spacing.")
