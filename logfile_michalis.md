LOGFILE for the project of COMPARATIVE GENOMICS

Author: LAZAR Michalis-Daniel



\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*  Thursday 27 November 2025 \*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*



I started by doing a brief literature review to understand the contenxt our project relating to past research.



Found the following platform that reunites all the possible resources on Danio rerio: https://zfin.org/



Also learned about \*Ohnologs\* which is a special subtype of paralog that arose from a whole-genome duplication event as opposed to

a small local/segmental duplication and named after Susumu Ohno who proposed that genome duplication 

was key in vertebrate evolution.



Found the following paper on patterns of ohnolog retention in 6 salmon species, maybe we could look into it for Danio rerio.

(https://doi.org/10.1002/ece3.9994)



I reproduce the TD with the file Erine provided to me (blast protein db all-against-all):



>head zebrafish\_blast\_results\_clean.txt

ENSDARP00000134140.2    ENSDARP00000134140.2    100.000 888     0       0       1       888     1       888     0.0     1846    888     888

ENSDARP00000134140.2    ENSDARP00000152482.1    99.105  782     7       0       1       782     1       782     0.0     1613    888     784

ENSDARP00000134140.2    ENSDARP00000127253.1    78.076  894     186     5       2       888     3       893     0.0     1459    888     893

ENSDARP00000134140.2    ENSDARP00000137017.1    72.306  863     225     5       26      882     28      882     0.0     1283    888     886

ENSDARP00000134140.2    ENSDARP00000147426.1    68.343  875     261     5       15      885     11      873     0.0     1246    888     874

ENSDARP00000134140.2    ENSDARP00000115823.2    68.343  875     261     5       15      885     11      873     0.0     1246    888     874

ENSDARP00000134140.2    ENSDARP00000150373.1    63.043  920     281     9       1       885     2       897     0.0     1153    888     898

ENSDARP00000134140.2    ENSDARP00000096257.3    31.835  267     161     7       607     869     144     393     1.11e-32        131     888     399

ENSDARP00000134140.2    ENSDARP00000146359.1    32.197  264     158     7       607     866     142     388     1.52e-32        130     888     396

ENSDARP00000134140.2    ENSDARP00000110495.2    32.197  264     158     7       607     866     142     388     1.52e-32        130     888     396



Reminder!The columns are the following:

1: qprot

2: sprot

3: pident (%identity)

4: alen (alignment length)

5: mismatches

6: gaps

7: qstart

8: qend

9: sstart

10: send

11: evalue 

12: bitscore

13: qlen (query length)

14: slen (subject length)



Coverage on each protein: cov\_q = alen/qlen, cov\_s = alen/slen



I need to decide on some concrete thresholds:



**High stringency (HS): very high-confidence paralogs**



pident ≥ 80

cov\_q ≥ 0.8 and cov\_s ≥ 0.8

evalue ≤ 1e-20



**Medium stringency (MS):**



pident ≥ 60

cov\_q ≥ 0.6 and cov\_s ≥ 0.6

evalue ≤ 1e-10



**Low stringency (LS): more explorative (twilight-zone of homology)**



pident ≥ 30

cov\_q ≥ 0.3 and cov\_s ≥ 0.3

evalue ≤ 1e-5

REMEMBER! I must remove self-hits in all three datasets


awk 'BEGIN{OFS="\\t"}

{

&nbsp;   qid = $1; sid = $2;

&nbsp;   if (qid == sid) next;              #to remove self-hits



&nbsp;   pid  = $3 + 0;                     #to awk coerce them into numbers instead of strings

&nbsp;   alen = $4 + 0;

&nbsp;   qlen = $13 + 0;

&nbsp;   slen = $14 + 0;

&nbsp;   eval = $11 + 0;



&nbsp;   cov\_q = alen / qlen;

&nbsp;   cov\_s = alen / slen;



&nbsp;   if (pid >= 80 \&\& cov\_q >= 0.8 \&\& cov\_s >= 0.8 \&\& eval <= 1e-20) #playing with the thresholds for each!

&nbsp;       print qid, sid, pid, alen, qlen, slen, eval, $12;

}' zebrafish\_blast\_results\_clean.txt > zebrafish\_HS.txt





>wc -l zebrafish\_\*.txt

&nbsp;   50783 zebrafish\_HS.txt

&nbsp; 2210443 zebrafish\_LS.txt

&nbsp;  204684 zebrafish\_MS.txt

&nbsp; 5201789 zebrafish\_blast\_results\_clean.txt



>head zebrafish\_HS.txt

ENSDARP00000134140.2    ENSDARP00000152482.1    99.105  782     888     784     0       1613

ENSDARP00000136350.1    ENSDARP00000152373.1    97.619  462     456     462     0       936

ENSDARP00000133086.1    ENSDARP00000131686.2    99.169  361     361     361     0       746

ENSDARP00000138273.1    ENSDARP00000015862.8    82.727  330     330     350     0       572

ENSDARP00000138224.1    ENSDARP00000135021.1    96.768  495     495     495     0       967

ENSDARP00000126495.2    ENSDARP00000059226.4    81.905  210     210     210     1.15e-130       366

ENSDARP00000152911.1    ENSDARP00000136147.1    93.312  942     952     934     0       1704

ENSDARP00000152911.1    ENSDARP00000145216.1    98.865  793     952     814     0       1596

ENSDARP00000152911.1    ENSDARP00000150944.1    95.244  799     952     855     0       1550

ENSDARP00000152911.1    ENSDARP00000156672.1    91.039  837     952     842     0       1535





Now we convert each filtered file into an undirected weighted graph file for clustering:

we will have node1 (protein ID A), node2 (protein ID B) and weight=bitscore. We need to collapse (A,B)

and (B,A) into a single entry.



awk 'BEGIN{OFS="\\t"}

{

&nbsp;   a = $1; b = $2; score = $8;    #bitscore

&nbsp;   if (a == b) next;

&nbsp;   if (a < b) print a, b, score; #lexiconographically compare and canonize order

&nbsp;   else       print b, a, score;

}' zebrafish\_HS.txt | sort -u > zebrafish\_HS\_forMCL.txt #deduplicate (collapse) duplicate lines



>head zebrafish\_HS\_forMCL.txt

ENSDARP00000000042.9    ENSDARP00000136868.2    1702

ENSDARP00000000042.9    ENSDARP00000137512.2    1708

ENSDARP00000000069.7    ENSDARP00000149136.1    724

ENSDARP00000000160.9    ENSDARP00000127358.1    768

ENSDARP00000000382.7    ENSDARP00000081489.4    4867

ENSDARP00000000382.7    ENSDARP00000081489.4    4895

ENSDARP00000000488.6    ENSDARP00000124141.2    303

ENSDARP00000000803.9    ENSDARP00000005224.6    1732

ENSDARP00000000803.9    ENSDARP00000005224.6    1757

ENSDARP00000000803.9    ENSDARP00000115698.1    1548





>wc -l zebrafish\_\*\_forMCL.txt

&nbsp;  35750 zebrafish\_HS\_forMCL.txt

&nbsp;1807006 zebrafish\_LS\_forMCL.txt

&nbsp; 164621 zebrafish\_MS\_forMCL.txt



I went on Galaxy Europe and used the following tool: Markov Cluster Algorithm for graphs (Galaxy Version 22.282+galaxy0)

First I used the default parameter for "Inflation" which determines the cluster granularity:2.0 (how would the families=clusters change if I change the parameter?)



MCL: “random walks stay inside dense groups” idea.



We have "one cluster per line" format in the output files (each line corresponds to one family of proteins).



To make work easier, we want a file with 2 columns (proteinID clusterID)



nl -ba #number all lines exactly as they appear

nl -ba MCL\_families\_HS.tabular | awk '{

&nbsp; fam = $1

&nbsp; for (i = 2; i <= NF; i++) {

&nbsp;   g = $i

&nbsp;   sub(/\[()]/, "", g) #Rm any trailing space if any 

&nbsp;   if (g != "") print g, fam

&nbsp; }

}' > zebrafish\_HS\_MCL\_families.txt





>head zebrafish\_HS\_MCL\_families.txt

ENSDARP00000105975.3 1

ENSDARP00000106009.2 1

ENSDARP00000111449.1 1

ENSDARP00000114164.2 1

ENSDARP00000115007.2 1

ENSDARP00000116114.3 1

ENSDARP00000117957.2 1

ENSDARP00000121725.3 1

ENSDARP00000125125.1 1

ENSDARP00000125164.1 1



Now I will do some histograms in Rstudio.



The histograms reveal that almost all families are tiny. Everything is piled up in the first bin.

Most "families" are just pairs: 



summary(as.numeric(sizes\_HS))

\# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 

\#  2.000   2.000   2.000   2.641   2.000 115.000 



3326 families of size 2



Max size = 115

Prob one ancestral gene has many very similar copies that still pass the strict criteria.

Might be interesting to look into this family!



It's getting late, will continue with more stats and exploration of the families in R later.

Can't wait to get to GO overrepresentation in Panther.



\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*  Friday 28 November 2025 \*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

Some more exploratory data analysis in R on my clusters (remember Inflation=2, might be interesting for someone, maybe Raul, to do the clusters for other values like 1.5)

In paralog analysis, the distribution of family sizes is what matters.

summary(as.numeric(sizes_HS))
max(as.numeric(sizes_HS)) #115
table(as.numeric(sizes_HS))[1:10]   # how many families of size 1,2,...,10
#
 Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  2.000   2.000   2.000   2.641   2.000 115.000 
[1] 115

   2    3    4    5    6    7    8    9   10 
3326  675  185   64   41   16   13    6   11 
  11 
   9 

#same for MS

Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  2.000   2.000   2.000   3.163   3.000 264.000 
[1] 264

   2    3    4    5    6    7    8    9   10 
3301 1068  433  206  121   61   44   33   23 
  11 
  11 
#SAME FOR LS
 Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  2.000   2.000   3.000   5.175   5.000 916.000 
[1] 916

   2    3    4    5    6    7    8    9   10 
1964  999  602  382  229  158  118   73   62 
  11 
  39 

Took me a while to review a bit how R works... made a lot of different figures.

For us, HS represents the most biologically credible paralog families, so I will focus on them.

In HS we expect to have preserved WGD ohnolog pairs, recent paralogs, well-conserved duplicates

MS and LS are mostly just for sensitivity and completeness checks.

Biological interpretation of the histogram for HS (raw linear):

zebrafish retains many 2-gene ohnolog pairs (3326) that had strong conservation in coding sequence (since they passed they strict filters imposed in the beginning of the pipeline)

limited expansion beyond pairs (well, there is some but the size is noticeably different)

Note: we observe the power-law behavior from the overall histogram (which shows up a lot in nature; patern behind extreme events): a lot of families of paralogs of size 2 and a long heavy- tail (right skewed)

I got a bit bored of this part (might be back) so I proceeded to do the orthologous enrichment.

I selected two datasets for HS (families with size 2 and size more than 10 in accordance with TD part II)
We have the following biological question: 
What functions are enriched among genes that live in big, highly redundant families as opposed to genes that only have a single very close paralog?

and so I went to https://pantherdb.org/ and chose "Statistical overrepresentation test", uploaded my file, selected the organism (Denio rerio), annotation set (PANTHER GO-Slim Biological Process) and lastly used the default  gene list for Danio rerio.

We are using Fisher’s Exact for each GO term and FDR correction for multiple testing correction.

review moment: for each GO term T, the test asks "Among the genes in my list, is T more common than expected if I had picked genes at random from the whole genome?"

If genes were random, the fraction with T in my list ≈ fraction with T in whole genome!

Null hypothesis: “Term T is not enriched; your list is just a random sample from the genome regarding T.”

Fisher computes the probability of seeing a contingency table as extreme or more extreme than (a,b,c,d) under that null (hypergeometric distribution)

exported the results as table

now let's do the same for familiies of size 10 or more.

tired, will analyse the results tomorrow.

I JUST SAW RAUL PUSHED THE FILES WITH THE OG OVERREPRESENTATION ANALYSIS BEFORE ME! at least I practiced a bit! but I wish we allocate tasks (and communicate) more efficiently so we can cover the breadth of the project in these high stake times!!!

