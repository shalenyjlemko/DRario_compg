
# LOGFILE for the project of COMPARATIVE GENOMICS

**Author**: Raul DURAN DE ALBA

**Student number**: 20245534

**Group members**: Erine BENOIST, Denys BURYI and Michalis-Daniel LAZAR

*******************************

## Before Thursday November 27th, 2025

Found FISHTEDB in class, proposed Danio rero as it’s a model organism used to compare in humans immunological issues. We downloaded the TEs in from the website. Denys started a GitHub. I added the templates to work on the article and presentation on Latex. We met once in Mike’s and I house in order to work on the project and defined further advancement.

### Problems

We thought we had to little TE´s. Only 2275 TE entries. We still wanted to continue using Danio rerio as it is a model organism.

*******************************

## Thursday November 27th, 2025

We met at Mike’s and I place to work on the project. We understood that the TE’s present on FishTEDB are referencing TE families (2279 TE families). There are around 4535283 TEs positions pairs in the reference Genome. Obtained by:

``` bash
wc -l /Users/raulduran/Downloads/Danio_rerio.age
```

We uploaded the documents to GitHub and faced commit issues. I had to rebase which complicated a previous push on the statistics of the reference genome Denys had committed previously. I will try to run the proteome database creation that Erine already executed or just continue with the blast she produced. I also runned RepeatMasker in the France Galaxy Server on the Danio rerio reference genome to obtain more TEs positions with the TE families obtained from FishTEDB. We will compare both datasets.

From the Blast produced by Erine, we assue the following columns in order to procede:

1. qseqid: query sequence id
2. sseqid: subject sequence id
3. pident: percentage of identical positions
4. length: alignment length (sequence overlap)
5. mismatch: number of mismatches
6. gapopen: number of gap openings
7. qstart: start of alignment in query
8. qend: end of alignment in query
9. sstart: start of alignment in subject
10. send: end of alignment in subject
11. evalue: expect value
12. bitscore: bit score
13. qlen: Query sequence length
14. slen: Subject sequence length

*******************************

## Friday November 28th, 2025

We continued working on the project. I created a stringent dataset from the blast results provided by Erine.

```bash
BLAST="zebrafish_blast_results_clean.txt"

awk -v OFS='\t' '!/^#/{
  ql=$13; sl=$14; qs=$7; qe=$8; ss=$9; se=$10;
  if(ql>0 && sl>0){
    qc=((qe>=qs?qe-qs+1:qs-qe+1)/ql)*100;
    sc=((se>=ss?se-ss+1:ss-se+1)/sl)*100;
    if(qc>=80 && sc>=80 && $3>=70 && $11<=1e-10 && $12>=80)
      printf "%s\t%.2f\t%.2f\n",$0,qc,sc
  }}' "$BLAST" > dataset_stringent.txt

```

Where we filter the blast results to have at least 80% of coverage in both query and subject, at least 70% of identity, e-value less than 1e-10 and bit score greater than 80.

WE verify the number of entries in the stringent dataset:

``` bash
wc -l dataset_stringent.txt
#   116178 dataset_stringent.txt
```

After that we can obtain a list to input the sets of genes for clustering into MCL. Where we habve both id's and the bitscore.

``` bash
awk -v OFS='\t' '{print $1,$2,$12}' dataset_stringent.txt > zebrafish_mcl_input.txt
```

After that we can run MCL in EU Galaxy server. We obtain the .tabular file with the clusters. We verify how many clusters we have by counting the columns of the file and the lines (members of each cluster):

``` bash
wc -l MCL_stringent.tabular
#   21358 MCL_stringent.tabular
```

``` bash
awk -F'\t' '{print NF}' MCL_stringent.tabular | sort -nu | tail -n 1
#   188
```

We have 188 clusters and the biggest cluster has 21358 members.
The commands were saved in the Blast_processing.bash script.
We will continue working on the analysis of the clusters by retrieving the list of genes for families size 2 and the families with more than 10 members.

First we retrieve the families of size 2 by extracting columns with size 2, then we will obtain a list with all of the rows inside those column in one column and then sort according to alphabetical order and remove duplicates:

``` bash
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
' MCL_stringent.tabular MCL_stringent.tabular > family_size2.tabular
```

Then we obtain the list of genes inside those families of size 2:

``` bash
tr '\t' '\n' < family_size2.tabular \
  | sed '/^$/d' \
  | sort -u \
  > genes_size2_unique.txt
```

Now we continue with the families of size greater than 10 members.

``` bash
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
' MCL_stringent.tabular MCL_stringent.tabular > family_size_ge10.tabular
```

Then we obtain the list of genes inside those families of size greater than 10:

``` bash
tr '\t' '\n' < family_size_ge10.tabular \
  | sed '/^$/d' \
  | sort -u \
  > genes_size_ge10_unique.txt
```

The commands were saved into the MCL_processing.bash script.

We perform a preliminary functional analysis of those genes using PantherDB with the lists obtained. The analysis is the PANTHER Overrepresentation Test (Released 20240807) with the following parameters for both lists:

- Annotation Data Set: PANTHER GO-Slim Biological Process
- Reference List: Danio rerio (all genes in database)
- Statistical Test: Fisher's Exact
- Correction: FDR

We will continue with the functional analysis of those genes on R later.

*******************************
