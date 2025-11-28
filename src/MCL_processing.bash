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
FNR==NR {
    # First pass: count non-empty cells in each column
    for (i = 1; i <= NF; i++) {
        if ($i != "") cnt[i]++
    }
    if (NF > maxNF) maxNF = NF
    next
}
!init {
    # Determine which columns to keep (size == 2)
    for (i = 1; i <= maxNF; i++) {
        if (cnt[i] == 2) keep[i] = 1
    }
    init = 1
}
{
    # Print only the selected columns (still tab-separated)
    first = 1
    for (i = 1; i <= NF; i++) {
        if (keep[i]) {
            if (!first) printf FS
            printf "%s", $i
            first = 0
        }
    }
    print ""
}
' "$MCL" "$MCL" > family_size2.tabular

tr '\t' '\n' < family_size2.tabular \
  | sed '/^$/d' \
  | sort -u \
  > genes_size2_unique.txt

wc -l genes_size2_unique.txt

# Now for family size >=10

awk -F'\t' '
FNR==NR {
    # First pass: count non-empty cells in each column
    for (i = 1; i <= NF; i++) {
        if ($i != "") cnt[i]++
    }
    if (NF > maxNF) maxNF = NF
    next
}
!init {
    # Determine which columns to keep (size >= 10)
    for (i = 1; i <= maxNF; i++) {
        if (cnt[i] > 10) keep[i] = 1
    }
    init = 1
}
{
    # Print only the selected columns (still tab-separated)
    first = 1
    for (i = 1; i <= NF; i++) {
        if (keep[i]) {
            if (!first) printf FS
            printf "%s", $i
            first = 0
        }
    }
    print ""
}
' "$MCL" "$MCL" > family_size_ge10.tabular

tr '\t' '\n' < family_size_ge10.tabular \
  | sed '/^$/d' \
  | sort -u \
  > genes_size_ge10_unique.txt

wc -l genes_size_ge10_unique.txt
