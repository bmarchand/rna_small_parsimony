import matplotlib.pylab as plt
import json
import os

FAMILIES = [fam.split('/')[-1].split('.')[0] for fam in os.listdir('resources/tree_files/')]

# Filtering families to retain only interesting ones, i.e:
#       1. For which there is agreement between the seed file and tree file compositions
#       2. With non empty structures at the leaves

with open('divergence.json') as f:
    d = json.load(f)

not_nice = []
for family in FAMILIES:
    max_num_bps = 0
    try:
        if len(open('results/small_parsimony_input_rfam/'+family+'_leaf_annotated.tab').readlines()) > 10000000000000000:
            not_nice.append(family)
            continue
        for line in open('results/small_parsimony_input_rfam/'+family+'_leaf_annotated.tab').readlines():
#            if line.find('UNKNOWN') >= 0:           
#                not_nice.append(family)
#                break       
            if line.find(':') >= 0:
                struct = line.split(' ')[-1].rstrip('\n')
                max_num_bps = max(max_num_bps, struct.count('('))
        if max_num_bps==0 or len(struct) > 100 or d[family]==0:
            not_nice.append(family)
    except FileNotFoundError:
        not_nice.append(family)

NICE_FAMILIES = [fam for fam in FAMILIES if fam not in not_nice]
print(len(NICE_FAMILIES),'nice families')
print('excluding',len([k for k,v in d.items() if v==0]),'with same str everywhere')

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
axs[0,0].scatter([d2[f] for f in NICE_FAMILIES],[d1[f] for f in NICE_FAMILIES],s=6, c='black')
axs[0,0].errorbar(np.mean([d2[f] for f in NICE_FAMILIES]),np.mean([d1[f] for f in NICE_FAMILIES]), 
                  xerr=np.std([d2[f] for f in NICE_FAMILIES]),
                  yerr=np.std([d1[f] for f in NICE_FAMILIES]),
                  fmt='s',c='red')
axs[0,0].set_xlabel('number of bps\n with IL_ILC')
axs[0,0].set_ylabel(r'number of bps''\n'r'with RF_$\emptyset$')
axs[0,0].plot([0,50],[0,50],color='r')
#axs[0,0].xaxis.set_label_position('top')

axs[0,1].scatter([d3[f] for f in NICE_FAMILIES],[d1[f] for f in NICE_FAMILIES],s=6, c='black')
axs[0,1].errorbar(np.mean([d3[f] for f in NICE_FAMILIES]),np.mean([d1[f] for f in NICE_FAMILIES]), 
                  xerr=np.std([d3[f] for f in NICE_FAMILIES]),
                  yerr=np.std([d1[f] for f in NICE_FAMILIES]),
                  fmt='s',c='red')
axs[0,1].plot([0,50],[0,50],color='r')
#axs[0,1].xaxis.set_label_position('top')

axs[1,1].scatter([d3[f] for f in NICE_FAMILIES],[d2[f] for f in NICE_FAMILIES],s=6, c='black')
axs[1,1].errorbar(np.mean([d3[f] for f in NICE_FAMILIES]),np.mean([d2[f] for f in NICE_FAMILIES]), 
                  xerr=np.std([d3[f] for f in NICE_FAMILIES]),
                  yerr=np.std([d2[f] for f in NICE_FAMILIES]),
                  fmt='s',c='red')
axs[1,1].plot([0,50],[0,50],color='r')
axs[1,1].set_ylabel('number of bps\n with IL_ILC')
axs[1,1].set_xlabel('number of bps\n with RF_ILC')

axs[0,2].scatter([d4[f] for f in NICE_FAMILIES],[d1[f] for f in NICE_FAMILIES],s=6, c='black')
axs[0,2].errorbar(np.mean([d4[f] for f in NICE_FAMILIES]),np.mean([d1[f] for f in NICE_FAMILIES]), 
                  xerr=np.std([d4[f] for f in NICE_FAMILIES]),
                  yerr=np.std([d1[f] for f in NICE_FAMILIES]),
                  fmt='s',c='red')
axs[0,2].plot([0,50],[0,50],color='r')

axs[1,2].scatter([d4[f] for f in NICE_FAMILIES],[d2[f] for f in NICE_FAMILIES],s=6, c='black')
axs[1,2].errorbar(np.mean([d4[f] for f in NICE_FAMILIES]),np.mean([d2[f] for f in NICE_FAMILIES]), 
                  xerr=np.std([d4[f] for f in NICE_FAMILIES]),
                  yerr=np.std([d2[f] for f in NICE_FAMILIES]),
                  fmt='s',c='red')
axs[1,2].plot([0,50],[0,50],color='r')

axs[2,2].scatter([d4[f] for f in NICE_FAMILIES],[d3[f] for f in NICE_FAMILIES],s=6, c='black')
axs[2,2].errorbar(np.mean([d4[f] for f in NICE_FAMILIES]),np.mean([d3[f] for f in NICE_FAMILIES]), 
                  xerr=np.std([d4[f] for f in NICE_FAMILIES]),
                  yerr=np.std([d3[f] for f in NICE_FAMILIES]),
                  fmt='s',c='red')
axs[2,2].plot([0,50],[0,50],color='r')
axs[2,2].set_ylabel('number of bps \n with RF_ILC')
axs[2,2].set_xlabel('number of bps 'r'with IL_$\emptyset$')

fig.savefig('figures/scater_plot_roots.pdf',bbox_inches='tight')
plt.show()
