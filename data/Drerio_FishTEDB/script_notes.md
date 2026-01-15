

**TE insertion time** was obtained using the formula T = K/2r
(35), which is widely used for the estimation of TE inser-
tion time (36, 37), wherein T represents the insertion time,
K represents the Kimura distance-based copy divergence of
TEs, and r represents the nucleic acid substitution rate. To
obtain the K value, we used the method (https://github.
com/4ureliek/Parsing-RepeatMasker-Outputs) developed by
Kapusta et al. (38) to convert the divergence value in the
RepeatMasker result (.out file) to the K value.

**HOW TO EXPLAIN THE TE INSERTION PEAK**
Within teleost fish, elevated TE activity has been shown to coincide with species radiations in salmonids: [sauce](https://link.springer.com/article/10.1186/1471-2164-8-422#Sec12)

*DNA transposons gather mutations over time, resulting in vertical inactivation [33]. Eventually all that remains are relics of the original sequence; the majority of DNA transposons in genomes are no more than fragments and altered sequences. Therefore, most instances are remnants, and only those that have **recently been active can be found as intact or near intact copies**.*

Speciation and Adaptive Radiation: TE expansions might drive or accompany speciation events by promoting genetic diversification. Active transposons can cause insertional mutations and chromosomal rearrangements that lead to reproductive isolation
link.springer.com
link.springer.com
. The zebrafish lineage may have undergone a rapid speciation phase (perhaps related to the diversification of cyprinid fishes in the Paleogene), during which TEs were simultaneously active. In support, bursts of TE activity coinciding with species radiations have been documented in other fish groups: e.g., explosive adaptive radiations of African cichlids correlate with elevated TE content and new insertions
link.springer.com
. Under this scenario, TEs could have contributed to genomic innovations in zebrafish (novel regulatory elements or genetic novelty) that facilitated speciation. Conversely, the genomic instability from active TEs might itself have fragmented populations (through chromosomal incompatibilities), spurring speciation. This is intertwined with the “genomic shock” idea since rapid environmental adaptation often underpins radiations.

**CHROMOSOME 4**
[reference](https://www.nature.com/articles/nature12111)
The long arm of chromosome 4 is unique among zebrafish genomic regions, owing to its relative lack of protein-coding genes and its extensive heterochromatin. Chromosome 4 is known to be late-replicating and hybridization studies suggest that genomic copies of 5S ribosomal DNA (rDNA), which are not notably present on any other chromosome, are scattered along the long arm at high redundancy18. Immediately after the presumed centromere at approximately 24 megabases (Mb), the sequence landscape (Fig. 1 and Supplementary Fig. A4) shows a remarkable increase in repeat content, which continues through to the telomere of the long arm. At approximately 27 Mb, the otherwise uniform presence of the satellite repeat SAT-2 on the long arm ends abruptly. This location is also the starting point of uniform MOSAT-2 distribution, a satellite repeat that is nearly absent from all other chromosomes but highly enriched on the long arm of chromosome 4. 


## 1. Version 1 – how they *built* the TE library and annotated genomes


### (1) Discover candidate TEs

Three complementary strategies, all run on each fish genome:

* **De novo** (no prior knowledge)

  * **RepeatModeler → RECON + RepeatScout** on the genome
  * Filter out tandem repeats (>25% TRF coverage)
  * BlastX vs Swiss-Prot – remove hits that are actually protein-coding genes
  * Remove tRNA / rRNA using tRNAscan-SE + Rfam (Blastn filters)

* **Structure-based** (for retrotransposons)

  * **LTR retrotransposons**: LTR_STRUC + MGEScan-LTR

    * Keep only elements with: 2 similar LTRs, internal region, correct flanking di/tri-nucleotides, TSDs.
  * **Non-LTR retrotransposons**: MGEScan-non-LTR (HMM-based).

* **Homology-based** (for DNA transposons)

  * **TESeeker** with a curated set of 257 fish transposase proteins from RepBase/NCBI.
  * It BLASTs, assembles hits (CAP3), aligns them (ClustalW2) and outputs consensus contigs.
  * They keep only the highest-quality consensus contigs.

So v1 *library building* = union of de novo + structure-based + homology-based consensus TEs.

### (2) Classify and de-duplicate consensus TEs

* **Classification**

  * Some tools (RepeatModeler, MGEScan, TESeeker) already assign superfamilies (hAT, Gypsy, etc.).
  * For “unknowns”:

    * **REPCLASS** → tries to assign them to **superfamilies**.
    * Remaining unclassified: **TEclass (SVM)** → at least assign TE **order** (DNA, LTR, LINE, SINE).

* **Redundancy removal (v1)**

  * They **merged all consensus sequences into one big set**, which is highly redundant.
  * Use **CD-HIT-EST (90% identity)** to collapse near-duplicates, **but**:

    * Nested TEs (DNA transposon inside an LTR) are a problem: if you cluster everything together, the shorter TE can be “eaten” by the longer one.
    * So they **de-redundantize per superfamily**: run CD-HIT separately for hAT, Gypsy, L2, etc., to avoid nested elements killing each other.
  * “Unknown” consensus sequences: align back to the genome (Blast; ≥85% identity, ≥50% coverage) and keep only if **copy number ≥ 3** (to remove spurious one-offs).

Result of v1 pipeline: for each species, a **clean, non-redundant, classified TE library** that is then used as the RepeatMasker library to annotate the genome.

---

## 2. Version 2 (the one you use) – what changed

Two key differences for your project:

### (A) How they treat redundancy in v2

* In v1 they were quite aggressive about collapsing redundancy in the consensus library (per superfamily).
* In v2 they **simplify that step**: they de-redundantize at the **coarser level** (DNA / LTR / LINE / SINE / Unknown) instead of per superfamily, AND:
* **Critically:** they **do not de-duplicate the RepeatMasker hits themselves** anymore.

  > “We have not removed the redundancy of the RepeatMasker result because it can retain more comprehensive data… This increases the probability of identifying TE–gene interactions.”

So:

* v1: Non-redundant TE library → non-redundant annotation.
* v2: Still a cleaned library, but **RepeatMasker output is kept “raw”** (all overlapping fragments, subcopies, etc.), so the database has **all TE fragments with full RepeatMasker statistics** (percdiv, percdel, percins, etc.).

That’s exactly the detailed GFF + age file you’re using.

### (B) TE insertion time estimation (v2 only)

New in v2: they put a **date** on each TE copy.

Pipeline:

1. **K (Kimura distance)**

   * Take RepeatMasker `.out` divergence (%),
   * Convert to Kimura 2-parameter **K** using the **Kapusta parseRM pipeline**.

2. **r (substitution rate) per lineage**

   * Whole-genome alignments with **LASTZ + UCSC axtChain/ChainNet + MULTIZ**, zebrafish as reference.
   * Use **PHAST**:

     * `msa_view` to extract 4-fold degenerate sites,
     * `phyloFit` to estimate a time-calibrated phylogenetic tree; branch lengths in substitutions per site.
   * For each lineage, compute **root-to-tip substitution rate** and divide by **622.6 Myr** (chordate–vertebrate divergence) → **r** (subs/site/Myr).
   * These r values are tabulated in their Supplementary Table S2.

3. **Compute insertion time**

   * For each TE copy:
     [
     T = \frac{K}{2r}
     ]
   * Store this as “Insertion time estimation = X Mya” in the GFF/age files.

So v2 = v1 library + **per-copy ages**, + **full, non-collapsed RepeatMasker annotation**.

---

## 3. How to show this on slides

### Slide 1 – “FishTEDB v1: how they built the TE library”

**Layout:**

* Title: *“FishTEDB v1 – TE discovery and library construction”*

* Left side: simple flowchart, three boxes feeding into one.

  * Box 1: **De novo TE discovery**

    * bullet: RepeatModeler (RECON + RepeatScout)
    * sub-bullets: filter tandem repeats (TRF); filter proteins (BlastX vs Swiss-Prot); filter tRNA/rRNA (tRNAscan-SE + Rfam).

  * Box 2: **Structure-based (retrotransposons)**

    * LTR_STRUC + MGEScan-LTR (keep only elements with 2 LTRs, IR, TSDs)
    * MGEScan-non-LTR for non-LTR retrotransposons.

  * Box 3: **Homology-based (DNA transposons)**

    * TESeeker with 257 fish transposase proteins.

  * Arrows from the three boxes into a big box: **“Union of candidate TE consensus sequences”**.

* Right side: bullet list titled “Classification & redundancy removal”:

  * REPCLASS → superfamilies; TEclass → DNA/LTR/LINE/SINE.
  * CD-HIT (90% identity), de-redundant **per superfamily** (handle nested TEs).
  * Keep only “Unknown” elements with copy number ≥ 3.

### Slide 2 – “FishTEDB v2: what changed (our data)”

Split into two horizontal panels:

**Top half: “Annotation strategy”**

* Bullet points:

  * *“TE library construction is essentially the same, but:”*
  * Redundancy removal simplified at the level DNA/LTR/LINE/SINE/Unknown instead of fine superfamilies.
  * **Crucial change:** *RepeatMasker output is no longer collapsed*.

    * “All individual hits + percdiv/percdel/percins are kept → better TE–gene relationship analysis.”

You can illustrate with two little cartoons:

* v1: TE hits overlapping, then a “merge” icon → one cleaned feature.
* v2: same overlapping hits but kept separate.

**Bottom half: “New in v2 – dating each TE copy”**

* Show the formula big in the middle:

  [
  T = \frac{K}{2r}
  ]

* Underneath, three arrows:

  * **K (Kimura distance)** – from parseRM on RepeatMasker `.out`.
  * **r (substitution rate)** – from multi-species alignments (LASTZ + MULTIZ + PHAST) calibrated with the phylogenetic tree and 622.6 Myr.
  * **T (Mya)** – stored as “Insertion time estimation” and used in your plots.

Add a tiny diagram of a tree with zebrafish highlighted and “r” on the branch.

---

If you present it like that, your professor sees:

* You understand **where the TE library comes from** (not a black box).
* You understand **why v2 is richer**: more detailed annotation + a proper evolutionary dating of each TE copy, which is exactly what you exploited in the TE age distributions and TE–gene analyses.
