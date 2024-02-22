Scripts used for the MEMORI longitudinal genetic analyses in EAC during neoadjuvant treatment


# CNA analyses: 

1_Bin_sequenza_segments_into_100kB_bins.R \
This script bins Sequenza files with allele-specific copy number profiles into 100kb long bins.

2_Plottingsegments_CNS.R \
This script plots copy number state for the entire genome per sample and annotates clinical information (patient ID, Responsiveness and Timepoint).

3_CNA_percentage.R \
This script plots the percentage of exome with CNS>3, CNS=3, >2:0, 2:0, 1:0, 1:1 in REP and NRP at each timepoint.

4_CNA_clonality.R \
This script plots the fragment size of clonal, subclonal and private CNAs in REP and NRP. 


# SNV analyses: 

5_Filtering_SNVs.R \
This script filters somatic SNVs, which passed Annova filtering criteria, show > 5 total reads in tumour and normal control and > 3 variant reads in tumour sample. 

6_Mobster.R \
This script computes the CCF using Mobster package. Based on their CCF SNVs are classified as clonal or subclonal.

7_Clonality_SNVs.R \
This script plots total, clonal and subclonal mutational burden in different clinical subgroups.

8_Phylogenetic_tree.R \
These scripts run phylogenetic analyses and plot phylogenetic trees for patients with 3 or more samples.

9_Calling_COSMIC_signatures.R \
This script runs SigDeconstruct to call COSMIC signatures.

10_Plotting_Cosmic_signatures.R \
This script plots COSMIC signatures proportions in clinical subgroups and proportions of base changes in treatment-naive, chemo-induced and radiochemo induced SNVs. 

11_dNdS_analysis.R \
This script runs dNdS analysis on naive, chemo-induced and radiochemo induced SNVs. 


# Intogen driver genes:
 
12_Plot_genetic_alterations_driver.R \
This script plots CNA, SNVs and indels in 108 Intogen driver genes for oesophageal and gastric cancer 

12_Plot_genetic_alterations_driver.R \
This script quantifies  Intogen driver genes with genetic alterations (CNAs, SNVs or indels).


# Neoantigens:

14_Filtering_neoantigenic_SNVs.R \
This script filters neoantigenic SNVs, which passed Annova filtering criteria, show > 5 total reads in tumour and normal control and > 3 variant reads in tumour sample. 

15_Clonality_neoantigenic_SNVs.R \
This script plots total, clonal and subclonal neoantigenic SNV burden in different clinical subgroups.

16_Neoantigenic_SNV_expression.R \
This script plots the mean expression of genes coding for neoantigens during treatment.

17_Correlation_expression_neoSNV_immune_cells.R
This script correlates T-cell infiltrates with neoantigen expression. 
