# Scripts for generating main and supplementary figures of Leventhal et al. in preparation

```0_nd_screen_pathway_enrichment.R```: Code for identifying enriched pathways in the neurodegeneration screen hits with gProfiler, table S1

```1_GTEX_age_tpm_regression.R ```: Code for analyzing the changes in neurodegeneration screen hits in GTEX, generates figures 2A-D and supplementary figures 1 and 5

```2_single_cell_preprocess_doublet_filter.py```: Code for preprocessing Mathys et al. 2019 single nucleus RNA-seq data for batch-correction, filtering and differential expression

```3_single_cell_enrichment_analysis.R```: Code for batch correcting, clustering and differential expression analysis of Mathys et al. 2019 single-nucleus RNA-seq data for figures 2E and S2 

```4_genome_track_prepare.R```: Generate ChIP-seq and motif tracks for visualization in figure 3B

```5_plot_ChIP_seq_motif_tracks.py```: Visualize ChIP-seq peaks and TF binding motifs that overlap with eQTLs of interest, as visualized in figure 3B 

```6_eGene_pathway_supplementary_figure.R```: Scripts for identifying pathway signatures in the temporal cortex pyramidal neuron data and generating figure S3
```7_piumet_preprocess.R```: Script for preparing metabolomic data for use in PiuMet

```8_drosophila_proteomics_phosphoproteomics.R```: Script for generating figures 3C-E and figure S4, writes a differential phosphoproteomics file for use in later scripts

```9_multi_omic_prize_generation.R```: Script for preparing input into OmicsIntegrator in order to generate network visualized in Figure 4

```10_process_RNA_salmon```: Generate RNA counts from raw files using Salmon

```11_NPC_DeSeq.R```: Code for differential expression and gene set enrichment analysis from neural progenitor cell RNA-seq after CRISPRi knockdown of genes of interest. Used to generate figures 5D, 5E, 6F, 6G, S6, S7, Tables S8 and S9

```12_generate_hold_out_prizes.py```: Generate prize files for omics integrator removing one prize type at a time

```13_run_oi_datasource_drop.py```: Wrapper script for running ```oi_run.py```, which launches OmicsIntegrator given a text file listing the prize files, in this case each missing one input type to assess its importance

```14_betweenness_analysis.R```: Code to calculate betweenness of neurodegeneration screen hits and nodes annotated for involvement in DNA damage-related processes

```15_hyperparameter_statistics.py```: Script to generate supplementary figure 1, which depicts the distributions of robustness and specificity of nodes before and after filtering.

```16_GNN_grid_search.py```: Code for performing a hyperparameter sweep for a graph neural network model

```17_GNN_run.py```: Run the Graph Neural Network using the hyperparameters inferred from ```16_GNN_grid_search.py```

```18_network_propagation_analysis```: Analysis for determining the effects of gene expression after gene knockdown in neural progenitor cells as a function of network distance

```19_PCNA_quantification_fig6h```: Code for generating figure 6h that quantifies PCNA foci after CSNK2A1 knockdown in Drosophila neurons

```20_figure4_circos_plot.R```: Code for plotting figure 4a, which represents the PCSF-derived network as a chord diagram

```21_cell_type_proportion_TCPY```: Code for generating Extended Data figure 3a, which shows the expression of Allen Brain Atlas marker genes of brain cell types in the temporal cortex pyramidal neuron RNA-seq data
