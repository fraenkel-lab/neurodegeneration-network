import pandas as pd, os, sys
import re

def main():
    prizefile=pd.read_excel('Supplementary\ table\ 7.xlsx',sheet=1)
    
    print(len(prizefile.index))
    print(prizefile['source'].unique())
    
    prizefile_copy=prizefile.copy()
    
    print(prizefile_copy['source'].unique())
    
    #rename data types
    prizefile_copy.loc[prizefile_copy['source'].isin(['drosophila ABeta proteomics',
                                                     'drosophila tau proteomics']),'source']='drosophila_proteomics'
    print(prizefile_copy['source'].unique())
    prizefile_copy.loc[prizefile_copy['source'].isin(['abeta drosophila phosphoproteomics',
                                                     'tau drosophila phosphoproteomics']),
                                                     'source']='drosophila_phosphoproteomics'
    prizefile_copy.loc[prizefile_copy['source'].isin(['abeta drosophila metabolomics',
                                                     'tau drosophila metabolomics']),'source']='drosophila_metabolomics'
    prizefile_copy.loc[prizefile_copy['source'].isin(['GWAS locus','accessible GWAS locus']),'source']='human_GWAS_locus'
    prizefile_copy.loc[prizefile_copy['source'].isin(['male-specific APOE4 lipidomics',
                                                     'female-specific APOE4 lipidomics']),'source']='human_lipidomics'
    
    prizelist=pd.DataFrame()
    
    #loop through each data type and drop one source
    for feature in prizefile_copy['source'].unique():
        prizefile_tmp=prizefile_copy[~(prizefile_copy['source'].isin([feature]))]
        
        feature=re.sub(' ','_',feature)
        print(feature)
        
        prizefile_tmp.to_csv(f'{feature}_omit_prize.txt',sep='\t',index=False)
        
        temp_prizelist=pd.DataFrame({'prefix':[f'{feature}_omit']},index=[0])
        print(temp_prizelist)
        
        prizelist=pd.concat([prizelist,temp_prizelist])
        
    print(prizelist)
    
    ## to run in 13_rev_run_oi_datasource_drop.py
    prizelist.to_csv('omit_feature_filelist.txt',sep='\t',index=False)
    
if __name__=='__main__':
    main()