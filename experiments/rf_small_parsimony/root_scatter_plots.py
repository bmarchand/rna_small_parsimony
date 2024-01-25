import matplotlib.pylab as plt
import json
import os

FILTERED_RFAM = [line.rstrip('\n') for line in open('FILTERED_RFAM').readlines()]
VERY_FILTERED_RFAM = [line.rstrip('\n') for line in open('VERY_FILTERED_RFAM').readlines()]

print(len(FILTERED_RFAM),'families')

DIR1 = 'results/small_parsimony_results_rfam_fitch_rf/'
DIR2 = 'results/small_parsimony_results_rfam_median_heuristic/'
DIR3 = 'results/small_parsimony_results_rfam_c2_rf_median_heuristic/'
DIR4 = 'results/small_parsimony_results_rfam_unconstrained_il_median_heuristic/'

def value_dict(DIR):
    d = {}
    for fname in os.listdir(DIR):
        family = fname.split('_')[0]

        root_struct = open(DIR+fname).readlines()[0].split(' ')[-1].rstrip('\n')

        d[family] = root_struct.count('(')

    return d


d1 = value_dict(DIR1)
d2 = value_dict(DIR2)
d3 = value_dict(DIR3)
d4 = value_dict(DIR4)

fig, axs = plt.subplots(3,3,sharex=True,sharey=True)

for i in range(1,3,1):
    for j in range(i):
        axs[i,j].axis('off')

import numpy as np
axs[0,0].scatter([d2[f] for f in FILTERED_RFAM],[d1[f] for f in FILTERED_RFAM],s=6, c='black')
axs[0,0].errorbar(np.mean([d2[f] for f in FILTERED_RFAM]),np.mean([d1[f] for f in FILTERED_RFAM]), 
                  xerr=np.std([d2[f] for f in FILTERED_RFAM]),
                  yerr=np.std([d1[f] for f in FILTERED_RFAM]),
                  fmt='s',c='red')
axs[0,0].set_xlabel('number of bps\n with IL_ILC')
axs[0,0].set_ylabel(r'number of bps''\n'r'with RF_$\emptyset$')
axs[0,0].plot([0,50],[0,50],color='r')
#axs[0,0].xaxis.set_label_position('top')

axs[0,1].scatter([d3[f] for f in FILTERED_RFAM],[d1[f] for f in FILTERED_RFAM],s=6, c='black')
axs[0,1].errorbar(np.mean([d3[f] for f in FILTERED_RFAM]),np.mean([d1[f] for f in FILTERED_RFAM]), 
                  xerr=np.std([d3[f] for f in FILTERED_RFAM]),
                  yerr=np.std([d1[f] for f in FILTERED_RFAM]),
                  fmt='s',c='red')
axs[0,1].plot([0,50],[0,50],color='r')
#axs[0,1].xaxis.set_label_position('top')

axs[1,1].scatter([d3[f] for f in FILTERED_RFAM],[d2[f] for f in FILTERED_RFAM],s=6, c='black')
axs[1,1].errorbar(np.mean([d3[f] for f in FILTERED_RFAM]),np.mean([d2[f] for f in FILTERED_RFAM]), 
                  xerr=np.std([d3[f] for f in FILTERED_RFAM]),
                  yerr=np.std([d2[f] for f in FILTERED_RFAM]),
                  fmt='s',c='red')
axs[1,1].plot([0,50],[0,50],color='r')
axs[1,1].set_ylabel('number of bps\n with IL_ILC')
axs[1,1].set_xlabel('number of bps\n with RF_ILC')

axs[0,2].scatter([d4[f] for f in FILTERED_RFAM],[d1[f] for f in FILTERED_RFAM],s=6, c='black')
axs[0,2].errorbar(np.mean([d4[f] for f in FILTERED_RFAM]),np.mean([d1[f] for f in FILTERED_RFAM]), 
                  xerr=np.std([d4[f] for f in FILTERED_RFAM]),
                  yerr=np.std([d1[f] for f in FILTERED_RFAM]),
                  fmt='s',c='red')
axs[0,2].plot([0,50],[0,50],color='r')

axs[1,2].scatter([d4[f] for f in FILTERED_RFAM],[d2[f] for f in FILTERED_RFAM],s=6, c='black')
axs[1,2].errorbar(np.mean([d4[f] for f in FILTERED_RFAM]),np.mean([d2[f] for f in FILTERED_RFAM]), 
                  xerr=np.std([d4[f] for f in FILTERED_RFAM]),
                  yerr=np.std([d2[f] for f in FILTERED_RFAM]),
                  fmt='s',c='red')
axs[1,2].plot([0,50],[0,50],color='r')

axs[2,2].scatter([d4[f] for f in FILTERED_RFAM],[d3[f] for f in FILTERED_RFAM],s=6, c='black')
axs[2,2].errorbar(np.mean([d4[f] for f in FILTERED_RFAM]),np.mean([d3[f] for f in FILTERED_RFAM]), 
                  xerr=np.std([d4[f] for f in FILTERED_RFAM]),
                  yerr=np.std([d3[f] for f in FILTERED_RFAM]),
                  fmt='s',c='red')
axs[2,2].plot([0,50],[0,50],color='r')
axs[2,2].set_ylabel('number of bps \n with RF_ILC')
axs[2,2].set_xlabel('number of bps 'r'with IL_$\emptyset$')

fig.savefig('figures/scater_plot_roots.pdf',bbox_inches='tight')
plt.show()
