import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

# Colors consistent with your flip levels
flip_levels = ["10%", "20%", "50%"]
legend_colors = ["#aec7e8", "#1f77b4", "#08306b"]

# Create colored patches
handles = [mpatches.Patch(color=c, label=lbl) for c, lbl in zip(legend_colors, flip_levels)]

# Make a clean standalone legend figure
fig, ax = plt.subplots(figsize=(2.8, 0.8))
ax.axis("off")

legend = ax.legend(
    handles=handles,
    title="Missing Rate",
    loc="center",
    ncol=3,
    fontsize=16,
    title_fontsize=18,
    frameon=False,
)

fig.patch.set_facecolor("white")
plt.tight_layout()
plt.savefig("figures_boxplots_final/legend_missing_rate.png", dpi=300, bbox_inches="tight", transparent=True)
plt.savefig("figures_boxplots_final/legend_missing_rate.svg", bbox_inches="tight", transparent=True)
plt.close(fig)

print("✅ Saved legend figure: legend_missing_rate.png / .svg")
