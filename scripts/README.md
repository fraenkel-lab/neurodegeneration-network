# Scripts for generating main and supplementary figures of Leventhal et al. in preparation

```0_nd_screen_pathway_enrichment.R```: Code for identifying enriched pathways in the neurodegeneration screen hits with gProfiler, table S1

```1_GTEX_age_tpm_regression.R ```: Code for analyzing the changes in neurodegeneration screen hits in GTEX, generates figures 2A-D and supplementary figures 1 and 5

```2_single_cell_preprocess_doublet_filter.py```: Code for preprocessing Mathys et al. 2019 single nucleus RNA-seq data for batch-correction, filtering and differential expression

```3_single_cell_enrichment_analysis.R```: Code for batch correcting, clustering and differential expression analysis of Mathys et al. 2019 single-nucleus RNA-seq data for figures 2E and S2 

```4_genome_track_prepare.R```: Generate ChIP-seq and motif tracks for visualization in figure 3B

```5_plot_ChIP_seq_motif_tracks.py```: Visualize ChIP-seq peaks and TF binding motifs that overlap with eQTLs of interest, as visualized in figure 3B 

```6_piumet_preprocess.R```: Script for preparing metabolomic data for use in PiuMet

```7_drosophila_proteomics_phosphoproteomics.R```: Script for generating figures 3C-E and figure S4, writes a differential phosphoproteomics file for use in later scripts

```8_multi_omic_prize_generation.R```: Script for preparing input into OmicsIntegrator in order to generate network visualized in Figure 4

```9_process_RNA_salmon```: Generate RNA counts from raw files using Salmon

```10_NPC_DeSeq.R```: Code for differential expression and gene set enrichment analysis from neural progenitor cell RNA-seq after CRISPRi knockdown of genes of interest. Used to generate figures 5D, 5E, 6F, 6G, S6, S7, Tables S8 and S9
