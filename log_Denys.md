**long time ago**
let up the repository

# Getting data

**27.11**
downloaded the genome of the same version as the published TE database. 
Analyzing the genome file
## installed software
installed samtools
#####

```
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
The following additional packages will be installed:
  libhts3t64 libhtscodecs2 libncurses6
Suggested packages:
  cwltool
The following NEW packages will be installed:
  libhts3t64 libhtscodecs2 libncurses6 samtools
0 upgraded, 4 newly installed, 0 to remove and 124 not upgraded.
Need to get 1240 kB of archives.
After this operation, 2999 kB of additional disk space will be used.
Do you want to continue? [Y/n] y
Get:1 http://archive.ubuntu.com/ubuntu noble/main amd64 libncurses6 amd64 6.4+20240113-1ubuntu2 [112 kB]
Get:2 http://archive.ubuntu.com/ubuntu noble/universe amd64 libhtscodecs2 amd64 1.6.0-1build1 [98.8 kB]
Get:3 http://archive.ubuntu.com/ubuntu noble/universe amd64 libhts3t64 amd64 1.19+ds-1.1build3 [437 kB]
Get:4 http://archive.ubuntu.com/ubuntu noble/universe amd64 samtools amd64 1.19.2-1build2 [593 kB]
Fetched 1240 kB in 0s (5663 kB/s)
Selecting previously unselected package libncurses6:amd64.
(Reading database ... 62243 files and directories currently installed.)
Preparing to unpack .../libncurses6_6.4+20240113-1ubuntu2_amd64.deb ...
Unpacking libncurses6:amd64 (6.4+20240113-1ubuntu2) ...
Selecting previously unselected package libhtscodecs2:amd64.
Preparing to unpack .../libhtscodecs2_1.6.0-1build1_amd64.deb ...
Unpacking libhtscodecs2:amd64 (1.6.0-1build1) ...
Selecting previously unselected package libhts3t64:amd64.
Preparing to unpack .../libhts3t64_1.19+ds-1.1build3_amd64.deb ...
Unpacking libhts3t64:amd64 (1.19+ds-1.1build3) ...
Selecting previously unselected package samtools.
Preparing to unpack .../samtools_1.19.2-1build2_amd64.deb ...
Unpacking samtools (1.19.2-1build2) ...
Setting up libhtscodecs2:amd64 (1.6.0-1build1) ...
Setting up libncurses6:amd64 (6.4+20240113-1ubuntu2) ...
Setting up libhts3t64:amd64 (1.19+ds-1.1build3) ...
Setting up samtools (1.19.2-1build2) ...
Processing triggers for man-db (2.12.0-4build2) ...
Processing triggers for libc-bin (2.39-0ubuntu8.6) ...
```
#### seqkit
```
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
The following NEW packages will be installed:
  seqkit
0 upgraded, 1 newly installed, 0 to remove and 124 not upgraded.
Need to get 7086 kB of archives.
After this operation, 16.9 MB of additional disk space will be used.
Get:1 http://archive.ubuntu.com/ubuntu noble-updates/universe amd64 seqkit amd64 2.3.1+ds-2ubuntu0.3 [7086 kB]
Fetched 7086 kB in 0s (19.1 MB/s)
Selecting previously unselected package seqkit.
(Reading database ... 62358 files and directories currently installed.)
Preparing to unpack .../seqkit_2.3.1+ds-2ubuntu0.3_amd64.deb ...
Unpacking seqkit (2.3.1+ds-2ubuntu0.3) ...
Setting up seqkit (2.3.1+ds-2ubuntu0.3) ...
Processing triggers for man-db (2.12.0-4build2) ...
```
### Exploratory data analysis

```
samtools faidx GCF_000002035.6_GRCz11_genomic.fna

cat GCF_000002035.6_GRCz11_genomic.fna.fai
```
>output
| Sequence name | Length   | Offset    | Bases per line | Bytes per line |
|---------------|----------|-----------|----------------|----------------|
| NC_007112.7   | 59578282 | 80        | 80             | 81             |
| NC_007113.7   | 59640629 | 60323171  | 80             | 81             |
| NC_007114.7   | 62628489 | 120709388 | 80             | 81             |
| NC_007115.7   | 78093715 | 184120814 | 80             | 81             |

```
grep '^NC_' GCF_000002035.6_GRCz11_genomic.fna.fai
```
| Sequence name | Length   | Offset     |
|---------------|----------|------------|
| NC_007112.7   | 59578282 | 80         |
| NC_007113.7   | 59640629 | 60323171   |
| NC_007114.7   | 62628489 | 120709388  |
| NC_007115.7   | 78093715 | 184120814  |
| NC_007116.7   | 72500376 | 263190781  |
| NC_007117.7   | 60270059 | 336597492  |
| NC_007118.7   | 74282399 | 397621007  |
| NC_007119.7   | 54304671 | 472832016  |
| NC_007120.7   | 56459846 | 527815576  |
| NC_007121.7   | 45420867 | 584981252  |
| NC_007122.7   | 45484837 | 630969961  |
| NC_007123.7   | 49182954 | 677023440  |
| NC_007124.7   | 52186027 | 726821262  |
| NC_007125.7   | 52660232 | 779659696  |
| NC_007126.7   | 48040578 | 832978262  |
| NC_007127.7   | 55266484 | 881619429  |
| NC_007128.7   | 53461100 | 937576826  |
| NC_007129.7   | 51023478 | 991706271  |
| NC_007130.7   | 48449771 | 1043367624 |
| NC_007131.7   | 55201332 | 1092423099 |
| NC_007132.7   | 45934066 | 1148314529 |
| NC_007133.7   | 39133080 | 1194822852 |
| NC_007134.7   | 46223584 | 1234445177 |
| NC_007135.7   | 42172926 | 1281246637 |
| NC_007136.7   | 37502051 | 1323946806 |
| NC_002333.2   | 16596    | 1700402753 |

```
seqkit stat GCF_000002035.6_GRCz11_genomic.fna
```

| File        | Format | Type | Num Seqs | Sum Length    | Min Length | Avg Length | Max Length |
|-------------|--------|------|----------|---------------|------------|------------|------------|
| our_file.fna| FASTA  | DNA  | 1,923    | 1,679,203,469 | 650        | 873,220.7  | 78,093,715 |

```
seqkit fx2tab -n -l GCF_000002035.6_GRCz11_genomic.fna | head

seqkit fx2tab -n -l GCF_000002035.6_GRCz11_genomic.fna > danRer11_all_lengths.tsv
```

```
seqkit fx2tab -n -l GCF_000002035.6_GRCz11_genomic.fna | grep '^NC_' > danRer11_chr_lengths.tsv
```
To calculate N50:
sort by length desc, accumulate, stop at half

```
total=$(awk -F'\t' '{sum += $2} END {print sum}' danRer11_all_lengths.tsv)

sort -t$'\t' -k2,2nr danRer11_all_lengths.tsv \
| awk -F'\t' -v total="$total" '
    BEGIN { half = total / 2 }
    {
        cum += $2
        if (cum >= half) {
            print "N50 = " $2 " bp, sequence: " $1
            exit
        }
    }'
```

N50 = 52186027 bp, sequence: NC_007124.7 Danio rerio strain Tuebingen chromosome 13, GRCz11 Primary Assembly


```
# Total length of all sequences
total_all=$(awk -F'\t' '{sum += $2} END {print sum}' danRer11_all_lengths.tsv)

# Total length of chromosomes only (NC_ sequences)
total_chr=$(awk -F'\t' '{sum += $2} END {print sum}' danRer11_chr_lengths.tsv)

# Calculate percentages
echo "Total genome length: $total_all bp"
echo "Chromosome length: $total_chr bp"
echo "Scaffold length: $((total_all - total_chr)) bp"
echo ""
awk -v chr="$total_chr" -v all="$total_all" 'BEGIN {
    chr_pct = (chr / all) * 100
    scaf_pct = 100 - chr_pct
    printf "Chromosomes: %.2f%%\n", chr_pct
    printf "Scaffolds: %.2f%%\n", scaf_pct
}'
```

>output
Total genome length: 1679203469 bp
Chromosome length: 1345118429 bp
Scaffold length: 334085040 bp

Chromosomes: 80.10%
Scaffolds: 19.90%

I want to perform the same analysis on the TEs, to do this:
```
awk 'BEGIN{OFS="\t"}
     {len = $5 - $4 + 1; print $1,$4,$5,len,$7,$9}' Danio_rerio.age \
  > danRer11_te_simple.tsv

```
Columns now:

1) chromosome (seqid)

2) start

3) end

4) length

5) strand

6) age_Myr


```
wc -l danRer11_te_simple.tsv
4535283

# Total TE bp
awk '{sum += $4} END {print sum}' danRer11_te_simple.tsv

# Per-chromosome TE count and total bp
awk '{cnt[$1]++; bp[$1]+= $4} 
     END{for (c in cnt) print c, cnt[c], bp[c]}' danRer11_te_simple.tsv \
  | sort
```

## TE analysis

**01.12**

Wrote and ran the script to extract TE data out of the author's sorted gff
script [here](/projects/compG/data/Drerio_FishTEDB/extract_tsv.py)

Now performing TE analysis

**02.12**

Resolved a bunch of conflicts on the repository, investigated the issue with LFS being full. 

Wrote per chromosome TE composition code inside the [notebook](data/Drerio_FishTEDB/perChr_TE.ipynb)

Wrote script to extract TE density as well as gene density on the chromosome body using sliding window method. Found an appropriate .gff3 [file](data/Drerio_FishTEDB/Danio_rerio.GRCz11.115.gff3) from ensembl to perform this. 

Generated graphs for all of this

**03.12**

Investigated the graphs from 02.12, performed per chromosome TE composition analysis. Chose the appropriate graphs to include in the presentation based on that.

**03.12 c.d.**
After the class we all discussed what we have been doing and what needs to be done. As a result we (me and Raul) are discussing the strategy about what data we have in the current duplicated genes we are analyzing, and what analyses can we do on it. Based on the distribution of the Ka/Ks ratios we agreed on splitting the gene families into 4 categories based on their Ka/Ks:
<= 0.1 -- likely housekeeping genes, under very high purifying selection
0.1 > x < 0.9 -- genes under the purifying selection (subfunctionalization)
0.9 => x < 1.1 -- genes under neutral selection
above 1.1 -- genes under adaptive selection

**04.12**
Did some research on the teleost fish for the presetntation, finished up the slides and figures. 
