import pandas as pd, os

def main():
    # dictionary to read files from GEO repository (number forthcoming)
	sample_ids=pd.read_table('auxiiliary_files/salmon_sample_list.txt')

	for idx,row in sample_ids.iterrows():
		os.system('salmon quant -i salmon_partial_sa_index -l A \
		-1 '+str(row['R1'])+'-2 '+ str(row['R2'])+ '-p 8 --validateMapping \
		salmon_perturbseq_outs/'+row['ID']+'_quant')
		
if __name__ == '__main__':
		main()
