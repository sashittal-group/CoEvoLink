import pandas as pd
from pathlib import Path
import cherryml

# adjust these
IN_CSV = "interaction_matrix_tree_leaves_beta_only.csv"
HOST_COL = "host_105273"       # <--- pick the host column you want to study
OUT_DIR = Path("cherry_msa_dir")
OUT_DIR.mkdir(exist_ok=True)

# load csv (first column is virus id)
df = pd.read_csv(IN_CSV, index_col=0)

if HOST_COL not in df.columns:
    raise SystemExit(f"{HOST_COL} not found in columns: try one of:\n{list(df.columns)[:10]} ...")

# series: virus -> 0/1
states = df[HOST_COL].astype(int)

# FASTA with 0/1 characters (not all tools accept these chars)
fasta_01 = OUT_DIR / f"{HOST_COL}_01.txt"
with open(fasta_01, "w") as fh:
    for virus, val in states.items():
        fh.write(f">{virus}\n{val}\n")

# FASTA recoded to two letters (A for 0, C for 1) — recommended for CherryML
fasta_AC = OUT_DIR / f"{HOST_COL}_AC.txt"
with open(fasta_AC, "w") as fh:
    for virus, val in states.items():
        letter = "C" if val == 1 else "A"
        fh.write(f">{virus}\n{letter}\n")

print("Wrote:", fasta_01, fasta_AC)
