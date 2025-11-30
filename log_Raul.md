
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

After that we can run MCL in EU Galaxy server. We used the settings Inflation 2.0,  We obtain the .tabular file with the clusters. We verify how many clusters we have by counting the lines of the file and the columns (members of each cluster):

``` bash
wc -l MCL_stringent.tabular
#   21358 MCL_stringent.tabular clusters
```

``` bash
awk -F'\t' '{print NF}' MCL_stringent.tabular | sort -nu | tail -n 1
#   188 members
```

We have 21358 clusters and the biggest cluster has 188 members.
The commands were saved in the Blast_processing.bash script.
We will continue working on the analysis of the clusters by retrieving the list of genes for families size 2 and the families with more than 10 members.

First we retrieve the families of size 2 by extracting lines with 2 columns, then we will obtain a list with all of the genes inside those lines in one column and then sort according to alphabetical order and remove duplicates:

``` bash
awk -F'\t' '
  {
    members = 0
    for (i = 1; i <= NF; i++) if ($i != "") members++
    if (members == 2) print
  }
' MCL_stringent.tabular > family_size2.tabular
```

Then we obtain the list of genes inside those families of size 2:

``` bash
tr '\t' '\n' < family_size2.tabular | sed '/^$/d' | sort -u > genes_size2_unique.txt
```

Now we continue with the families of size greater than 10 members.

``` bash
awk -F'\t' '
  {
    members = 0
    for (i = 1; i <= NF; i++) if ($i != "") members++
    if (members > 10) print
  }
' MCL_stringent.tabular > family_gt10.tabular

```

Then we obtain the list of genes inside those families of size greater than 10:

``` bash
tr '\t' '\n' < family_gt10.tabular | sed '/^$/d' | sort -u > genes_gt10_unique.txt
```

The commands were saved into the MCL_processing.bash script.

We perform a preliminary functional analysis of those genes using PantherDB with the lists obtained. The analysis is the PANTHER Overrepresentation Test (Released 20240807) with the following parameters for both lists:

- Annotation Data Set: PANTHER GO-Slim Biological Process
- Reference List: Danio rerio (all genes in database)
- Statistical Test: Fisher's Exact
- Correction: FDR

We will continue with the functional analysis of those genes on R later.

*******************************

## Saturday November 29th, 2025

I realized my mistake at processing the MCL output file. I was counting columns as the clusters. I fixed the scripts accordingly and re-ran the analysis. Modify my bash scripts accordingly.
Now I started working in how to obtain the Ka/Ks ratios for the gene families. I started by creating a list of pairs from the stringent blast results to input into ParaAT. I extracted the first two columns of the zebrafish_mcl_input.txt file:

```bash
cut -f1,2 zebrafish_mcl_input.txt > pairs_for_paraat.txt
```

Then I prepared the CDS extracted from [Ensembl](https://ftp.ensembl.org/pub/release-115/fasta/danio_rerio/cds/). I cleaned the fasta file to have only the sequence identifiers in the headers.

```bash
awk '/^>/{sub(/^>/,">", $1); print $1; next} !/^>/' Danio_rerio.GRCz11.cds.all.fa > drerio.clean.fasta
```

 I also prepered my conda environment with the required tools: EMBOSS (for transeq),MAFFT and KaKsCalculator2.

```bash
conda install bioconda::emboss
conda install bioconda::mafft
conda install kakscalculator2
```

Then I translated the CDS fasta into a protein fasta using transeq from EMBOSS:

```bash
transeq -sequence drerio.clean.fasta -outseq zebrafish.pep.fasta -table 1
```

Now I needed a tool to automate the KaKs calculation. I downloaded [ParaAT2.0](https://ngdc.cncb.ac.cn/tools/paraat) from China National Center for Bioinformation and unzipped it. I added them to src folder. I made the scripts executable and added the current directory to the PATH.

```bash
chmod +x ParaAT.pl
chmod +x Epal2nal.pl
export PATH="$(pwd):$PATH"
```

Finally I ran ParaAT with the following command:

```bash
perl ../../src/ParaAT2.0/ParaAT.pl \
  -h pairs_for_paraat.txt \
  -n drerio.clean.fasta \
  -a zebrafish.pep.fasta \
  -m mafft \
  -p processors.txt \
  -f axt \
  -k \
  -o ParaAT_output
```

I realized that I needed to create a processors.txt file with a single line containing the number 4 to indicate the number of processors to use. After I run it i realized a major issue. The paths shared contain only the Protein IDs and not the CDS IDs. I will need to create a mapping file between both IDs to be able to run ParaAT correctly.
I eliminated the version information as the GTF didn't had it in their IDs and then created the mapping file by downloading the Danio rerio GTF file from Ensembl and extracting the relevant information with the following command:

```bash
awk '{
  id1 = $1;
  id2 = $2;
  sub(/\.[0-9]+$/, "", id1);   # remove .2, .7, etc.
  sub(/\.[0-9]+$/, "", id2);
  print id1 "\t" id2;
}' pairs_for_paraat.txt > pairs_for_paraat.noversion.txt
```

Then I extracted the relevant information from the GTF file:

```bash
perl -ne '
  next if /^#/;                          # skip comments
  next unless /\tCDS\t/ && /protein_id/; # only CDS lines that have a protein_id
  my ($tx)   = /transcript_id "([^"]+)"/;
  my ($prot) = /protein_id "([^"]+)"/;
  print "$prot\t$tx\n" if $tx && $prot;
' Danio_rerio.GRCz11.115.gtf > prot2tx.map
```

Then I proceded to run the mapping between both files to obtain the CDS IDs for each Protein ID in the pairs file:

``` bash
awk 'NR==FNR{
        map[$1]=$2;        # map[protein] = transcript   (from prot2tx.map)
        next
     }
     {
        if (1 in map && $2 in map && map[$1] != map[$2])
            print map[$1] "\t" map[$2];
     }' prot2tx.map pairs_for_paraat.noversion.txt > pairs_tx.txt
```

Now I have to verify I removed the versioning information correctly and then I will be able to run ParaAT again with the correct pairs file.

``` bash
# pep: remove `_1`
awk '/^>/{sub(/^>/,"",$1); sub(/_.*/,"",$1); print ">"$1; next} !/^>/' \
    zebrafish.pep.fasta > zebrafish.tx.pep.fasta

# cds: ensure header is just first token
awk '/^>/{sub(/^>/,"",$1); print ">"$1; next} !/^>/' \
    Danio_rerio.GRCz11.cds.all.fa > drerio.tx.cds.fasta
```

I will run ParaAT now with the corrected files.

```bash
perl ../../src/ParaAT2.0/ParaAT.pl \
  -h pairs_tx.txt \
  -n drerio.tx.cds.fasta \
  -a zebrafish.tx.pep.fasta \
  -m mafft \
  -p processors.txt \
  -f axt \
  -k \
  -o ParaAT_output
```

It started at Sat Nov 29 20:37:58 2025. It has encountered multiple inconsistences between the pep and nuc sequences files (ERROR: inconsistency between the following pep and nuc seqs. Run bl2seq (-p tblastn) or GeneWise to see the inconsistency.). KaKs calculator started at Sat Nov 29 21:56:12 2025. I will have to check the output files later. Decided due to slow process to run ParaAT on the local machine of my teamate Michalis-Daniel Lazar who has a more powerful computer.

I find it outstanding that RandomMasker in Galaxy France is still running after so many hours. It started more than 24h ago. I will try to use Galazy EU server to run it again.

I managed to scrape the data from FishTEDB using a wget command to obtain the full data of TEs in Danio rerio including sequences. I added the data to our GitHub repository in the data/Drerio_FishTEDB folder.

*******************************

## Sunday November 30th, 2025
