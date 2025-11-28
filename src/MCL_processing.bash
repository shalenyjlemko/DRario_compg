#!/bin/bash

# ========================================================
# MCL clustering processing for Zebrafish Blast results
# ========================================================

MCL="MCL_stringent.tabular"

# Count number of lines in MCL file which corresponds the number of maximum members in the biggest family
wc -l "$MCL"

# Determine maximum number of columns in MCL file, which indicates the number of families
awk -F'\t' '{print NF}' "$MCL" | sort -nu | tail -n 1


# We obtain the gene lists for family size ==2 

awk -F'\t' '
  {
    members = 0
    for (i = 1; i <= NF; i++) if ($i != "") members++
    if (members == 2) print
  }
' MCL_stringent.tabular > family_size2.tabular

# Now the list of unique genes in families of size 2

tr '\t' '\n' < family_size2.tabular | sed '/^$/d' | sort -u > genes_size2_unique.txt

# We obtain the gene lists for family size >10

awk -F'\t' '
  {
    members = 0
    for (i = 1; i <= NF; i++) if ($i != "") members++
    if (members > 10) print
  }
' MCL_stringent.tabular > family_gt10.tabular

# Now the list of unique genes in families of size >10

tr '\t' '\n' < family_gt10.tabular | sed '/^$/d' | sort -u > genes_gt10_unique.txt