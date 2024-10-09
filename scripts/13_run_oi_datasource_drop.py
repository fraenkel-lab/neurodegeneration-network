import pandas as pd, os, sys

def main():
    ## input file should be omit_feature_filelist.txt generated in 12_generate_hold_out_prizes.py
    filelist = sys.argv[1]
    
    prefixes = pd.read_table(filelist)
    for i,row in prefixes.iterrows():
        ## oi_run.py is a utility script that runs OmicsIntegrator2's implementation of the
        ## Prize-Collecting Steiner Forest (PCSF)
        call = "python3 oi_run.py str(row["prefix"])+"_prize.txt "+str(row["prefix"])
        print(call)
        os.system(call)
        
if __name__ == '__main__':
    main()
