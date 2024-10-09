import pickle, numpy as np, pandas as pd, networkx as nx, warnings
warnings.filterwarnings('ignore')

#import matplotlib, matplotlib.pyplot as plt
#matplotlib.rcParams['font.sans-serif'] = "Arial"
#matplotlib.rcParams['font.family'] = "Arial"

#import seaborn as sns,
import OmicsIntegrator as oi
import sys

def main():
    infile=str(sys.argv[1])
    prizefile=str(sys.argv[2])
    outname=str(sys.argv[3])

    Ws = [1,3,6]
    Bs = [2,5,10]
    Gs = [3,5,6]

    params = {
        "noise": 0.1, 
        "dummy_mode": "terminals", 
        "exclude_terminals": False, 
        "seed": 1
    }

    #interactome
    graph = oi.Graph(infile, params)
    graph.prepare_prizes(prizefile)

    node_attributes_df = graph.node_attributes

    robustness_reps = 100
    specificity_reps = 100

    results = graph.grid_randomization(prizefile, Ws, Bs, Gs, robustness_reps, specificity_reps)

    with open(str(outname)+"_iref_randomization_results.pkl", "wb") as f: 
        pickle.dump(results, f)


if __name__ == '__main__':
        main()

