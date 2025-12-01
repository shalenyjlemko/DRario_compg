import gzip
import csv
from pathlib import Path

gff_path = Path("Danio_rerio.fasta.out.gff.sorted.gz")
out_path = Path("Danio_rerio_TE_table.tsv")

with gzip.open(gff_path, "rt") as fin, open(out_path, "w", newline="") as fout:
    w = csv.writer(fout, delimiter="\t")
    # header
    w.writerow([
        "chrom", "start", "end", "length", "strand",
        "TE_name", "TE_class", "PercDiv", "PercDel", "PercIns", "age_Myr"
    ])

    for line in fin:
        if line.startswith("#") or not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue

        chrom, source, feature, start, end, score, strand, phase, attrs = fields
        start_i = int(start)
        end_i = int(end)
        length = end_i - start_i + 1

        attr_dict = {}
        for item in attrs.split(";"):
            item = item.strip()
            if not item:
                continue
            if "=" not in item:
                continue
            key, val = item.split("=", 1)
            attr_dict[key] = val

        te_name = attr_dict.get("Target", "NA")
        te_class = attr_dict.get("Class", "NA")

        perc_div = attr_dict.get("PercDiv", "NA")
        perc_del = attr_dict.get("PercDel", "NA")
        perc_ins = attr_dict.get("PercIns", "NA")

        age_str = attr_dict.get("Insertion time estimation", "")
        age_myr = "NA"
        if age_str:
            # "34.75 Mya" -> 34.75
            age_myr = age_str.strip().split()[0]
            # basic sanity: must look like a float
            try:
                float(age_myr)
            except ValueError:
                age_myr = "NA"

        w.writerow([
            chrom, start_i, end_i, length, strand,
            te_name, te_class, perc_div, perc_del, perc_ins, age_myr
        ])

print("Wrote:", out_path)