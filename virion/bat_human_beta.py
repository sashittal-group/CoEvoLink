import pandas as pd

# --- Input files ---
matrix_file = "flipped_matrix_fn_hosts.csv"
host_tax_file = "filtered_host_taxonomy.csv"
virus_tax_file = "filtered_virus_taxonomy.csv"

# --- Read all files ---
mat = pd.read_csv(matrix_file, dtype=str)
host_tax = pd.read_csv(host_tax_file, dtype=str).fillna("")
virus_tax = pd.read_csv(virus_tax_file, dtype=str).fillna("")

# --- List of viruses to analyze ---
target_viruses = ["3418604", "3433757", "694009"]

# Clean up matrix to have virus IDs as base numbers (remove 'virus_' prefix)
mat["VirusTaxID"] = mat.iloc[:, 0].apply(lambda x: x.split("_", 1)[1] if "_" in x else x)

# Melt the matrix to long form: each row is (virus, host, value)
long_df = mat.melt(id_vars=["VirusTaxID"], var_name="HostTaxID", value_name="Value")

# Clean up host IDs (remove 'host_' prefix)
long_df["HostTaxID"] = long_df["HostTaxID"].apply(lambda x: x.split("_", 1)[1] if "_" in x else x)

# Keep only 1's and target viruses
long_df = long_df[(long_df["Value"] == "1") & (long_df["VirusTaxID"].isin(target_viruses))]

# --- Merge with taxonomy data ---
merged = (
    long_df.merge(host_tax, on="HostTaxID", how="left", suffixes=("", "_host"))
           .merge(virus_tax, on="VirusTaxID", how="left", suffixes=("_host", "_virus"))
)

# Optional: select useful columns only
cols = [
    "VirusTaxID", "Virus", "VirusGenus", "VirusFamily", "VirusOrder", "VirusClass",
    "HostTaxID", "Host", "HostGenus", "HostFamily", "HostOrder", "HostClass"
]
merged = merged[cols]

# --- Save result ---
merged.to_csv("bat_human_virus_host_links.csv", index=False)

print(f"Saved {len(merged)} virus-host links for {len(target_viruses)} target viruses to bat_human_virus_host_links.csv")
