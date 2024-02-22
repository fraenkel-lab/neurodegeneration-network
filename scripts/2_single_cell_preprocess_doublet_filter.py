import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import seaborn as sns
import scrublet as scr

def main():
    adata=sc.read('raw_ROSMAP/matrix.mtx.gz',cache=True)
    adata.X = adata.X.astype('float64')
    counts_matrix = adata.X.tocsc()
    genes_scrub = np.array(gene_ids[0].values)
    
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)

    scrub.plot_histogram()
    
    doublet_df=adata.obs[predicted_doublets]
    doublet_df.to_csv('doublet_cells.txt',sep='\t',index=False)
    
if __name__=='__main__':
    main()
