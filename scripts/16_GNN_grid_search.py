from torch_geometric.data import Data, DataLoader
from torch_geometric.nn import GCNConv, Set2Set, GATConv, Linear
from torch_geometric.explain import Explainer,GNNExplainer
import torch_geometric.transforms as T
import torch
import torch.nn.functional as F
import os
from tqdm import tqdm, trange

import matplotlib.pyplot as plt
import pandas as pd, numpy as np
from sklearn.model_selection import train_test_split
#from joblib import Parallel, delayed

## GCN class
class Net(torch.nn.Module):
    def __init__(self, num_features, dim=16, num_hid=8, num_classes=1):
        super(Net, self).__init__()
        self.conv1 = GCNConv(num_features, dim)
        self.conv2 = GCNConv(dim, num_hid)
        self.linear = Linear(num_hid,num_classes)

    def forward(self, x, edge_index, data=None):
        x = F.relu(self.conv1(x, edge_index))
        x = F.dropout(x, training=self.training)
        x = F.relu(self.conv2(x, edge_index))
        x = self.linear(x)
        return x

#map node to index
def get_idx_dir(df):
    mapping = {index: i for i, index in enumerate(df.index.unique())}

    return mapping

def to_onehot_matrix(categorical_series):
    #written in gemini
    """
    Converts a pandas Series of categorical values to a one-hot encoded feature matrix.
    
    Args:
      categorical_series (pd.Series): A pandas Series containing categorical values.
    
    Returns:
      torch.Tensor: A one-hot encoded feature matrix suitable for the GAT network.
    """
    # Get unique categories
    unique_categories = categorical_series.unique()
    num_categories = len(unique_categories)
    
    # One-hot encode using category codes
    category_codes = categorical_series.factorize(sort=True)[0]
    onehot_matrix = torch.zeros(len(categorical_series), num_categories,dtype=torch.float32)
    onehot_matrix.scatter_(1, torch.tensor(category_codes).unsqueeze(dim=1), 1,reduce='add')
    
    return onehot_matrix

def test_reg(model,data,mode='test'):
    model.eval()

    if mode=='test':
        log_logits=model(data.x,data.edge_index,data)
        test_loss=F.mse_loss(log_logits[data.test_mask],data.y[data.test_mask])
    elif mode=='validate':
        log_logits=model(data.x[data.val_mask],data.edge_index[data.val_mask],
                         data)
        test_loss=F.mse_loss(log_logits,data.y[data.val_mask])

    return test_loss

def run_model(in_graph,device,dim,epochs=100,lr=0.01,weight_decay=5e-3):
    model = Net(num_features=in_graph.num_features, dim=dim, num_classes=in_graph.num_classes).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

    t = trange(epochs, desc="Stats: ", position=0)
    
    for epoch in t:
    
        model.train()
    
        in_graph = in_graph.to(device)
        optimizer.zero_grad()
        log_logits = model(in_graph.x, in_graph.edge_index, in_graph)
    
        # Since the data is a single huge graph, training on the training set is done by masking the nodes that are not in the training set.
        #loss = F.nll_loss(log_logits[data.train_mask], data.y[data.train_mask])
        loss = F.mse_loss(log_logits[in_graph.train_mask], in_graph.y[in_graph.train_mask])
        loss.backward()
        optimizer.step()
    
        # validate
        #train_acc, test_acc = test(model, data)
        train_loss = loss
        test_loss = test_reg(model,in_graph,mode='test')
        
        #t.set_description('[Train_loss:{:.6f} Train_acc: {:.4f}, Test_acc: {:.4f}]'.format(loss, train_acc, test_acc))
        t.set_description('[Train_loss:{:.6f}, Test_loss {:.4f}'.format(train_loss,test_loss))

    results_df=pd.DataFrame({'train_loss_MSE':train_loss.item(),'test_loss_MSE':test_loss.item(),
                            'n_epochs':epochs,'learning_rate':lr,'weight_decay':weight_decay},index=[0])

    results_df.to_csv(
        f'GCN_loss_epochs_{epochs}_learning_rate_{lr}_weight_decay_{weight_decay}.txt',
                      sep='\t',index=False)

def main():
    ## format data for py geometric
    prize_file=pd.read_excel('Supplementary\ table\ 7.xlsx',sheet=1)
    interactome=pd.read_table('auxiliary_files/iref17_iref14_braInMap_ptau_metab_lipid.txt')

    nodes=pd.DataFrame({'id':np.unique(interactome[['protein1','protein2']].values)},
                   index=np.unique(interactome[['protein1','protein2']].values))

    node_df=nodes.merge(prize_file[['id','prize_val','source']],on='id',how='left')
    node_df['prize_val']=node_df['prize_val'].fillna(0)

    node_df.index=node_df['id']
    node_mapping=get_idx_dir(nodes)

    prot_list1=[node_mapping[v] for v in interactome['protein1']]
    prot_list2=[node_mapping[v] for v in interactome['protein2']]
    
    mapped_interacs=pd.DataFrame({'p1':prot_list1,'p2':prot_list2})
    
    edge_index=torch.tensor(mapped_interacs.to_numpy().T)

    ## format prize file so that the feature is the data source, dummified, and index is node name
    node_df['source']=node_df['source'].fillna('None')
    node_df.index=node_df['source']
    
    mapping=get_idx_dir(node_df)
    dtype_list=[mapping[v] for v in node_df['source']]

    ## convert x to one-hot matrix
    x=to_onehot_matrix(node_df['source'])

    edge_attr=torch.tensor(interactome['cost'].values,dtype=torch.float32).reshape(len(interactome.index),1)
    y=torch.tensor(node_df['prize_val'].values,dtype=torch.float32).reshape(len(node_df.index),1)
    y_grp=torch.tensor(pd.cut(node_df['prize_val'],bins=10,labels=False)).reshape(len(node_df.index),1)

    # make torch dat
    in_graph=Data(x=x,edge_attr=edge_attr,edge_index=edge_index,y=y,num_classes=1)
    in_graph.y_grouped=y_grp
    in_graph.train_mask,in_graph.test_mask=train_test_split(np.arange(len(y_grp)),test_size=0.2,
                                           shuffle=True,stratify=y_grp)

    epoch_list = [25,50,100]
    lr_list = [1e-2,1e-3,1e-4]
    weight_list = [1e-3,5e-3,1e-4,5e-4]
    dim = 16
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Parallel(n_jobs=4)(delayed(run_model)(in_graph,device,dim,epochs,
    #                                       lr,weight_decay) for epochs,lr,weight_decay in zip(epoch_list,lr_list,weight_list))

    for epochs in epoch_list:
        for lr in lr_list:
            for weight_decay in weight_list:
                run_model(in_graph,device,dim,epochs,lr,weight_decay)
                

if __name__=='__main__':
    main()

    