import networkx as nx
import pandas as pd, numpy as np
from axial import axial
import pickle
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import matplotlib,matplotlib.pyplot as plt
import seaborn as sns
import community as community_louvain
from networkx.algorithms import community
from glob import glob
import leidenalg as la
import igraph
mpl.rcParams['pdf.fonttype'] = 42

def summary(results,mode):
    
    steiner_series=[]
    prize_series=[]
    
    for paramstring, graphs in results.items():
        df=pd.DataFrame.from_dict(dict(graphs['augmented_forest'].nodes(data=True)),orient='index')

        if len(df.index)>0:
            steiner=df[df['terminal']]
            prize=df[~df['terminal']]

            steiner=steiner[mode].rename(paramstring)
            prize=prize[mode].rename(paramstring)

            steiner_series.append(steiner)
            prize_series.append(prize)
        
        
    steiner_node_summary_df = pd.concat(steiner_series, axis=1).fillna(0)
    prize_node_summary_df = pd.concat(prize_series, axis=1).fillna(0)
    
    return steiner_node_summary_df, prize_node_summary_df

def summary_single(results,mode):
    
    series=[]
    
    for paramstring, graphs in results.items():
        df=pd.DataFrame.from_dict(dict(graphs['augmented_forest'].nodes(data=True)),orient='index')
        
        if len(df.index)>1:
            if mode=='membership':
                df['membership']=1

            temp_df=df[mode].rename(paramstring)

            series.append(temp_df)
        
        
    node_summary_df = pd.concat(series, axis=1).fillna(0)
    
    return node_summary_df

def summary_membership_filter(results):
    
    series=[]
    
    for paramstring, graphs in results.items():
        df=pd.DataFrame.from_dict(dict(graphs['augmented_forest'].nodes(data=True)),orient='index')
        if len(df.index)>0:

            df.loc[(df['terminal'])|((df['robustness']>0.4)&(df['specificity']<0.4)),'filter_membership']=1

            temp_df=df['filter_membership'].rename(paramstring)

            series.append(temp_df)       
        
    node_summary_df = pd.concat(series, axis=1).fillna(0)
    
    return node_summary_df

def summary_terminal(results):
    
    series=[]
    
    for paramstring, graphs in results.items():
        df=pd.DataFrame.from_dict(dict(graphs['augmented_forest'].nodes(data=True)),orient='index')
        if len(df.index)>0:
        
            df['n_terminal']=0
            df.loc[df['terminal'],'n_terminal']=1

            temp_df=df['n_terminal'].rename(paramstring)

            series.append(temp_df)
        
        
    node_summary_df = pd.concat(series, axis=1).fillna(0)
    
    return node_summary_df

def gene_fishing(results,gene_list):
    
    series=[]
    
    for paramstring, graphs in results.items():
        df=pd.DataFrame.from_dict(dict(graphs['augmented_forest'].nodes(data=True)),orient='index')

        print(paramstring,df[df.index.isin(gene_list)].index)
        
        temp_df=pd.DataFrame()
        
        series.append(temp_df)
        
        
    node_summary_df = pd.concat(series, axis=1).fillna(0)
    
    return node_summary_df

def summary_epi(results,epi_list,prot_list):
    series=[]
    
    for paramstring, graphs in results.items():
        df=pd.DataFrame.from_dict(dict(graphs['augmented_forest'].nodes(data=True)),orient='index')
        
        temp_cluster=[]
        if len(df.index)>0:
            for cluster in df['louvain_clusters'].unique():
                if len(
                    df[(df['terminal'])&(df['louvain_clusters'].isin([cluster]))].index)>0:
                    temp_cluster.append(len(df[(df['source'].isin(prot_list))&(
                        df['louvain_clusters'].isin([cluster]))].index)/len(
                        df[(df['terminal'])&(df['louvain_clusters'].isin([cluster]))].index))
                else:
                    temp_df=pd.DataFrame({str(paramstring):0},index=[0])

            temp_df=pd.DataFrame({str(paramstring):np.mean(temp_cluster)},index=[0])

            series.append(temp_df)
        else:
            temp_df=pd.DataFrame({str(paramstring):0},index=[0])
            
    node_summary_df = pd.concat(series, axis=1).fillna(0)
    
    return node_summary_df

def pcsf_pareto(hyperparameter_dataframe,list_of_features):
    #input is a dataframe of PCSF evaluation metrics and a list of features to evaluate
    #make sure features are ordered with "cost" in mind, meaning that we want to prioritize
    #smaller values and de-prioritize large values (so instead of robustness we want 1-robustness, for example) 
    #returns Pareto efficient row subset of pts
    # sort points by decreasing sum of coordinates
    # code is largely ripped from this stackoverflow post: https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python
    pts=hyperparameter_dataframe[list_of_features].values
    pts = pts[pts.sum(1).argsort()[::-1]]
    # initialize a boolean mask for undominated points
    # to avoid creating copies each iteration
    undominated = np.ones(pts.shape[0], dtype=bool)
    for i in range(pts.shape[0]):
        # process each point in turn
        n = pts.shape[0]
        if i >= n:
            break
        # find all points not dominated by i
        # since points are sorted by coordinate sum
        # i cannot dominate any points in 1,...,i-1
        undominated[i+1:n] = (pts[i+1:] >= pts[i]).any(1) 
        # keep points undominated so far
        pts = pts[undominated[:n]]
    return pts
    
def return_optimal_hyperparams(hyperparameter_dataframe,list_of_features,pcsf_pareto_out):
    ## get which rows of the original hyperparameter dataframe (hyperparameter_dataframe)
    ## correspond to the optimal solutions from pcsf_pareto (pcsf_pareto_out)
    ## list_of_features should be the subset of features that were the inputs of pcsf_pareto
    ## returns filtered dataframe indicating which hyperparameters
    ## are optimal
    mask=[]
    for i,row in enumerate(pcsf_pareto_out):
        idx=np.where(np.equal(hyperparameter_dataframe[
            list_of_features].values,pcsf_pareto_out[i,:]))[0][0]
        
        mask.append(idx)
        
    return hyperparameter_dataframe[hyperparameter_dataframe.index.isin(mask)]

def jaccard_index_reference(reference_network,test_network_dict):
    ## determine similarity of network with that of a reference network
    ## reference_network: network to determine similarity
    ## test_network: results with networks to compare to reference
    
    X=nx.read_graphml(reference_network)
    
    series=[]
    for paramstring, graphs in test_network_dict.items():
        df=pd.DataFrame.from_dict(dict(graphs['augmented_forest'].nodes(data=True)),orient='index')
        
        if len(df.index)>0:
            Y=graphs['augmented_forest']

            intersec=nx.intersection(X,Y)

            jaccard_idx=len(intersec.nodes)/len(set(list(X.nodes()) + list(Y.nodes())))

            temp_df=pd.DataFrame({str(paramstring):[jaccard_idx]},index=[0])

            series.append(temp_df)
        
    res_df = pd.concat(series, axis=1).fillna(0)
    
    return res_df

def louvain_clustering(nxgraph,res,partition=None):

    nx.set_node_attributes(nxgraph, {node: {'louvain_clusters':str(cluster)} for node,cluster in community_louvain.best_partition(
        nxgraph,partition=partition,resolution=res).items()})
    
def main():
    ## make plots for Supplementary Figure 1 showing the distribution of robustness and specificity across predicted
    ## nodes
    
    ## parameter we used in the final network
    selected_hyperparameter='W_3.00_B_5.00_G_5.00'
    
    with open("oi_outs.pkl", "rb") as f:
        dat=pickle.load(f)
    
    robustness_df=summary_single(dat,'robustness')
    specificity_df=summary_single(dat,'specificity')
    
    robustness_df['id']=robustness_df.index
    specificity_df['id']=specificity_df.index

    robust_melt=pd.melt(robustness_df, id_vars=['id'], 
                        value_vars=robustness_df.columns[~(robustness_df.columns.isin(['id']))])
    spec_melt=pd.melt(specificity_df, id_vars=['id'], 
                        value_vars=specificity_df.columns[~(specificity_df.columns.isin(['id']))])

    plot_df=robust_melt.merge(spec_melt,on=['id','variable'])
    
    in_prizes=pd.read_excel('Supplementary\ table\ 7.xlsx',sheet=1)
    
    ## get robustness and specificity of the predicted nodes
    robustness_steiner=robustness_df[~(robustness_df.index.isin(in_prizes['id']))]
    specificity_steiner=specificity_df[~(specificity_df.index.isin(in_prizes['id']))]
    
    ## make df for plotting 
    robustness_steiner['id']=robustness_steiner.index
    specificity_steiner['id']=specificity_steiner.index

    robust_melt=pd.melt(robustness_steiner, id_vars=['id'], 
                        value_vars=robustness_steiner.columns[~(robustness_steiner.columns.isin(['id']))])
    spec_melt=pd.melt(specificity_steiner, id_vars=['id'], 
                        value_vars=specificity_steiner.columns[~(specificity_steiner.columns.isin(['id']))])

    plot_df=robust_melt.merge(spec_melt,on=['id','variable'])
    print(plot_df.head())
    
    font = {'family' : 'sans serif',
        'weight' : 'normal',
        'size'   : 20}
    plt.rc('font', **font)
    plt.rcParams['figure.figsize'] = [12, 10]

    g=sns.histplot(plot_df[plot_df['variable'].isin([selected_hyperparamter])], x="value_x",col_wrap=9)
    plt.ylabel('Number of nodes')
    plt.xlabel('Robustness')
    plt.savefig('robustness_histogram_pre_filter.pdf',
               bbox_inches='tight')
    
    plt.gcf().clear()
    
    g=sns.histplot(plot_df[plot_df['variable'].isin([selected_hyperparamter])], x="value_y",col_wrap=9)
    plt.ylabel('Number of nodes')
    plt.xlabel('Specificity')
    plt.savefig('specificity_histogram_pre_filter.pdf',
               bbox_inches='tight')
    
    plt.gcf().clear()
                           
    ## filter on robustness, specificity, and membership in a connected component
    plot_count_df_filt=plot_df[(plot_df['id'].isin(in_prizes['id']))|
                                ((plot_df['value_x']>0.4)&(plot_df['value_y']<0.4))]

    ## limit to nodes with degree >0
    connected_df=pd.DataFrame()
    for paramstring, graphs in dat.items():
        nodes = [node for (node, val) in graphs['augmented_forest'].subgraph(
            plot_count_df_filt.loc[plot_count_df_filt['variable'].isin([paramstring]),
                                   'id'].values).degree() if val>0]

        temp_df=pd.DataFrame({'id':nodes,'hyperparameter':paramstring})

        connected_df=pd.concat([connected_df,temp_df])    
    
    ## limit to selected hypeparameter and plot distribution of robustness and specificity
    connected_df_merge=connected_df[connected_df['hyperparameter'].isin(
                           ['W_3.00_B_5.00_G_5.00'])].merge(in_prizes,on='id',how='left')
    connected_df_merge_nodup=connected_df_merge.drop_duplicates(subset=['id'])
    print(connected_df_merge_nodup)
    
    font = {'family' : 'sans serif',
        'weight' : 'normal',
        'size'   : 28}
    plt.rc('font', **font)
    plt.rcParams['figure.figsize'] = [12, 10]

    connected_df_rob=connected_df_merge_nodup.merge(plot_df,on='id').drop_duplicates(subset=['id'])

    sns.histplot(connected_df_rob,x='value_x',bins=50)
    plt.xlabel('Robustness')
    plt.ylabel('Number of nodes')
    plt.savefig('robustness_in_network.pdf',bbox_inches='tight')
    
    
    plt.gcf().clear()
    
    font = {'family' : 'sans serif',
        'weight' : 'normal',
        'size'   : 28}
    plt.rc('font', **font)
    plt.rcParams['figure.figsize'] = [12, 10]

    connected_df_rob=connected_df_merge_nodup.merge(plot_df,on='id').drop_duplicates(subset=['id'])

    sns.histplot(connected_df_rob,x='value_y',bins=50)
    plt.xlabel('Specificity')
    plt.ylabel('Number of nodes')
    plt.savefig('specificity_in_network.pdf',bbox_inches='tight')
    
    ## calculate modularity of network with Leiden and Louvain clustering
    
    # filtered network
    G=igraph.read('spec_abeta_tau_lipid_W_3.00_B_5.00_G_5.00_anno.graphml',format='graphml')
    partition=la.find_partition(G,la.ModularityVertexPartition)

    node_df = pd.DataFrame({attr: G.vs[attr] for attr in G.vertex_attributes()})
    ## annotate dataframe with leiden cluster
    node_df['leiden_clusters']=0
    for i,part in enumerate(partition):
        node_df.loc[node_df.index.isin(part),'leiden_clusters']=i

    print(node_df.head())
    
    ## write to table, also calculate modularities (Q-values) of leiden and louvain algorithms
    node_df.to_csv('leiden_clustered_nodes.txt',sep='\t',index=False)

    #louvain and leiden quality value
    print(f'louvain modularity = {h.modularity([int(clu) for clu in node_df.louvain_clusters])}')
    print(f'leiden modularity = {h.modularity(node_df.leiden_clusters)}')
    
if __name__=='__main__':
    main()

    
    
    