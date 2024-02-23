from Bio import SeqIO
import gzip
import pandas as pd, numpy as np
import re
import numpy as np
import pandas as pd
import glob
import os
import scipy.stats as stats

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches
from matplotlib.ticker import FormatStrFormatter
from cycler import cycler

import qtl.annotation
import qtl.genotype as gt
import qtl.map
import qtl.plot
import qtl.locusplot

import style

def main():

    ##plot overlap of ChIP-seq peaks with identified transcription factor binding motifs
    #output from 4_genome_track_prepare.R
    tf_df = pd.read_csv('eqt_motif_search/motif_chip_seq_peaks_HLA_DRB.txt', sep='\t')
    motif_s = tf_df[tf_df['type']=='motif']#.iloc[0]
    tf_df = tf_df[tf_df['type']=='peak']
    
    tf_df = tf_df.sort_values('tf_factor')

    ax = qtl.plot.setup_figure(4, 1.25, xspace=[1, 0.75], yspace=[0.5, 0.25])
    h = 0.66
    for _,r in tf_df.iterrows():
        d = r['stop'] - r['start'] + 1
        ax.add_patch(patches.Rectangle((r['start'], r['tf_factor']-1-h/2), d, h, 
        facecolor=[0.75]*3, zorder=-10))

    for i,motis in motif_s.iterrows():
        r = motis
        d = r['stop'] - r['start'] + 1
        ax.add_patch(patches.Rectangle((r['start'], r['tf_factor']-1-h/2*0.66), 
        d, h*0.66, facecolor=[0.2]*3, zorder=-10))
        ax.set_yticks(np.arange(tf_df.shape[0]))
        ax.set_yticklabels(tf_df['tf'])
        ylim = [-1+h/2, tf_df.shape[0]-h/2]

    ax.plot([32579035, 32579035], ylim, '--', c="#714A9D")
    ax.set_ylim(ylim)
    b = 0
    ax.set_xlim([tf_df['start'].min()-b, tf_df['stop'].max()+b])
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=4))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    if 0:
        ax.xaxis.tick_top()
        ax.spines['bottom'].set_visible(False)
    else:
        ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', which='both', direction='out', labelsize=10, length=0)
    ax.tick_params(axis='x', labelsize=10)
    ax.set_xlabel('Position (b)', fontsize=14)
    ax.set_ylabel('ChIP-seq peaks', fontsize=12)

    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['pdf.fonttype'] = 'truetype'
    plt.rcParams['ps.fonttype'] = 'truetype'
    plt.savefig('Figure_3B_bottom_tracks.pdf',width=14,height=8.67)
    
if __name__=='__main__':
    main()