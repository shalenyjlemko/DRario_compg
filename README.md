# DRario_compg
A project for indentification and annotation of duplicated genes and transposable elements in Danio Rerio

## Authors

    Erine BENOIST and Michaelis-Daniel LAZAR - Data Processing & BLAST

    Raul DURAN DE ALBA and Denys BURYI - TE Analysis


## Duplicated genes 

The part on the Duplicated Genes is mainly based on the history of duplication of Danio rerio.
An Teleost Genome Duplication (TGD/3R) event happened a round 320-350 millions years ago and the impact is that zebrafish have more genes (~26,000) than humans (~20,000) largely due to this event.
It can be proved by one of the paper (presence of 7 Hox clusters (vs 4 in humans)).

Our Project's Goal: To distinguish between these ancient "Ohnologs" (via low-stringency BLAST) and recent tandem duplications (via high-stringency BLAST).

### Erine parts on this 
Data Preprocessing & BLAST Results

1. Data Source & Cleaning I started by downloading the official Danio rerio proteome (GRCz11) from Ensembl. The raw data contained every known isoform for every gene, which would have created "fake" duplicates in our analysis. To fix this, I wrote a Python script to filter the file, retaining only the longest isoform for each gene to create a non-redundant canonical dataset.

2. Removing Bias (Mitochondria) I also filtered out mitochondrial proteins from the dataset. Since mitochondrial DNA evolves at a different rate and has a different inheritance pattern than the nuclear genome, keeping these entries would have skewed our whole-genome statistics.

3. Quality Control (QC) Before running the alignment, I generated histograms to check the data quality. The protein length distribution looks healthy (peaking ~400aa). Crucially, the pipeline successfully handled "giant" proteins like Titin (>29,000aa), confirming that we haven't lost large data points during the filtration.

4. All-vs-All BLASTp I have finished running the All-vs-All BLASTp (aligning the ~26,000 filtered proteins against themselves). I used specific formatting parameters to include query and subject lengths in the output, which is required for the coverage calculations we need to do next.

5. Accessing the Data The final BLAST results file (zebrafish_blast_results.txt.gz) is uploaded to our GitHub repository.

    Note: It is compressed to save space. You will need to run gunzip on it before using it for the TE or Clustering analysis.

# GUYS PLEASE ADD WHAT YOU DID SO I CAN UNDERSTAND AND CONCLUDE

## Transposable elements

## References:

- The genome paper 
Howe, K., Clark, M. D., Torroja, C. F., Torrance, J., Berthelot, C., Muffato, M., Collins, J. E., Humphray, S., McLaren, K., Matthews, L., McLaren, S., Sealy, I., Caccamo, M., Churcher, C., Scott, C., Barrett, J. C., Koch, R., Rauch, G., White, S., . . . Barker, G. (2013). The zebrafish reference genome sequence and its relationship to the human genome. Nature, 496(7446), 498–503. https://doi.org/10.1038/nature12111

In this paper, they confirm that around 70% of human genes have at least one obvious zebrafish ortholog, and many human genes have two zebrafish counterparts (co-orthologs) due to the duplication event.

- The teleost genome duplication event 
Amores, A., Force, A., Yan, Y., Joly, L., Amemiya, C., Fritz, A., Ho, R. K., Langeland, J., Prince, V., Wang, Y., Westerfield, M., Ekker, M., & Postlethwait, J. H. (1998). Zebrafish hox Clusters and Vertebrate Genome Evolution. Science, 282(5394), 1711–1714. https://doi.org/10.1126/science.282.5394.1711

This paper proved the duplication by showing that while mammals have 4 Hox gene clusters, zebrafish have 7 (resulting from the duplication of the original 4, followed by the loss of one cluster).

