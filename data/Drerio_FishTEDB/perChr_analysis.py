import pandas as pd

te = pd.read_csv(
    "Danio_rerio_TE_table.tsv",
    sep="\t"
)

# to drop familes with NA age
# te = te[te["age_Myr"].notna() & (te["age_Myr"] != "NA")]

# clean age to float just in case
te["age_Myr"] = pd.to_numeric(te["age_Myr"], errors="coerce")

# split TE_class into broad type and superfamily
# e.g. "DNA/Harbinger" -> type="DNA", family="Harbinger"
te[["type", "family"]] = te["TE_class"].str.split("/", 1, expand=True)

# per-chr TE family composition (use bp, not just counts)
te["bp"] = te["length"].astype(int)

per_chr_family = (
    te.groupby(["chrom", "family"])["bp"]
      .sum()
      .reset_index()
)

per_chr_family.head()
