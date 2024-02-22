Longitudinal RNA expression analyses in MEMORI cohort

# Expression analyses: 

1_Normalisation_gene_expression.R \
This script filters samples with >1.5 million reads and normalises the raw gene count into count per million (cpm) expression.

2_Unsupervised_Hallmark_cancer_enrichment.R \
Enrichment analyses for Hallmark of cancer pathways \
Unsupervised clustering based on Hallmark of cancer gene sets \
Loading plots showing enrichment of hallmark of cancer pathways

3_Supervised_Hallmark_cancer_enrichment.R \
This scripts performs supervised hierarchical clustering based on hallmark of cancer pathway expression.

4_KEGG_GO_analysis.R \
This script performs KEGG and GO enrichment analyses between different clinical subgroups.

5_Plotting_KEGG_enrichment.R \
This script creates dot plots showing significant KEGG pathway enrichments between different clinical subgroups.

6_Driver_gene_expression.R \
This script plots the expression of Intogen drivers during treatment.

7_Cibersort.R \
This script plots immune cell types deconvoluted via CIBERSORT in different clinical subgroups.

8_ConsensusTME.R \
This script deconvolutes different types of immune cells using ConsensusTME and shows immune cell types in different clinical subgroups.
