import matplotlib.pyplot as plt
import numpy as np

# ==============================
# Data setup
# ==============================

datasets = ["cow", "gut", "water"]
methods = ["PHIST", "PBLKS", "WISH"]

# Their accuracies (1-hit)
their_1 = {
    "cow": [0.61, 0.015, 0.576],
    "gut": [0.43, 0.0119, 0.381],
    "water": [0.45, 0.0099, 0.554],
}

# Our accuracies (1-hit)
our_1 = {
    "cow": [0.61, 0.06, 0.576],
    "gut": [0.43, 0.024, 0.381],
    "water": [0.495, 0.03, 0.554],
}

# Their accuracies (2/3-hit)
their_23 = {
    "cow": [0.71, 0.106, 0.652],
    "gut": [0.63, 0.13, 0.50],
    "water": [0.61, 0.139, 0.752],
}

# Our accuracies (2/3-hit)
our_23 = {
    "cow": [0.73, 0.106, 0.71],
    "gut": [0.63, 0.19, 0.56],
    "water": [0.61, 0.178, 0.79],
}

# ==============================
# Plot styling
# ==============================

# Two distinct color schemes: one for 1-hit, one for 2/3-hit
colors = {
    "their_1": "#8da0cb",  # blue tone
    "our_1": "#4daf4a",    # green tone
    "their_23": "#fc8d62", # orange tone
    "our_23": "#e41a1c",   # red tone
}

# ==============================
# Plotting loop
# ==============================

for dataset in datasets:
    x = np.arange(len(methods))  # positions for PHIST, PBLKS, WISH
    bar_width = 0.18

    fig, ax = plt.subplots(figsize=(8, 5))

    # offsets for 1-hit and 2/3-hit groups
    offset_1 = -bar_width * 1.2
    offset_23 = bar_width * 1.2

    # --- 1-hit bars ---
    ax.bar(x + offset_1 - bar_width/2, their_1[dataset], width=bar_width,
           color=colors["their_1"], label="Other Method (1-hit)")
    ax.bar(x + offset_1 + bar_width/2, our_1[dataset], width=bar_width,
           color=colors["our_1"], label="CoEvoLink (1-hit)")

    # --- 2/3-hit bars ---
    ax.bar(x + offset_23 - bar_width/2, their_23[dataset], width=bar_width,
           color=colors["their_23"], label="Other Method (2/3-hit)")
    ax.bar(x + offset_23 + bar_width/2, our_23[dataset], width=bar_width,
           color=colors["our_23"], label="CoEvoLink (2/3-hit)")

    # Axis formatting
    ax.set_xticks(x)
    ax.set_xticklabels(methods, fontsize=11)
    ax.set_ylabel("Accuracy", fontsize=11)
    ax.set_ylim(0, 1)
    ax.set_title(f"{dataset.upper()} dataset", fontsize=14, fontweight="bold")
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    # Legend at upper right
    ax.legend(loc="upper right", fontsize=9, frameon=True)

    plt.tight_layout()

    # Save the figure
    plt.savefig(f"{dataset}_accuracy_comparison.svg", dpi=300)
    plt.savefig(f"{dataset}_accuracy_comparison.png", dpi=300)
    plt.close()
