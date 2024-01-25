import os
from rnadist.utils import tree_tab_file_parse

FILTERED_RFAM = [line.rstrip('\n') for line in open('FILTERED_RFAM').readlines()]
VERY_FILTERED_RFAM = [line.rstrip('\n') for line in open('VERY_FILTERED_RFAM').readlines()]

import json
with open('divergence.json') as f:
    divergence = json.load(f)

HIGHEST_DIVERGENCE = list(sorted(divergence.keys(), key=lambda x: -divergence[x]))[:60]

height_method = 'max'
metric='num_bps'
def height(node):
    if node.isleaf:
        if height_method=='num_desc':
            return 1
        else:
            return 0
    if height_method=='max':
        return max([height(child) for child in node.children])+1
    if height_method=='min':
        return min([height(child) for child in node.children])+1
    if height_method=='num_desc':
        return sum([height(child) for child in node.children])

def values_per_height_dict(DIR,family):

    for filename in os.listdir(DIR):
        print(filename)
        if  filename.split('_')[0]==family:
            lines = open(DIR+filename).readlines()
    
            phylo_T = tree_tab_file_parse(lines)
    
            queue = [phylo_T]

            d = {}
    
            while len(queue) > 0:
                node = queue.pop()
                h = height(node)
                if metric=='num_bps':
                    num = node.annotation.count('(')
                if metric=='num_loops':
                    num = num_loops(node.annotation)
                assert(node.annotation.count('(')==node.annotation.count(')'))
                try:
                    d[h].append(num)
                except KeyError:
                    d[h] = [num]
    
                for child in node.children:
                    queue.append(child)


    return [float(max(d[h])) for h in sorted(d.keys())]

DIR = 'results/small_parsimony_results_rfam_fitch_rf/'
DIR2 = 'results/small_parsimony_results_rfam_median_heuristic/'
DIR3 = 'results/small_parsimony_results_rfam_c2_rf_median_heuristic/'
DIR4 = 'results/small_parsimony_results_rfam_unconstrained_il_median_heuristic/'
DIR5 = 'results/small_parsimony_results_rfam_re_distance_heuristic/'

import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
c = sns.color_palette("light:b", as_cmap=True)

def plot(DIR,ax=None,cbar=False):
    list_dict = {}

    m = 0
    for family in HIGHEST_DIVERGENCE:
        list_dict[family] = values_per_height_dict(DIR,family)
        list_dict[family] = [v/max(list_dict[family]) for v in list_dict[family]]
        m = max(m, len(list_dict[family]))

    T = np.zeros((len(HIGHEST_DIVERGENCE),m))

    for i, family in enumerate(HIGHEST_DIVERGENCE):
        for j, v in enumerate(list_dict[family]):
            T[i,j] = v

#    qm = sns.heatmap(T,yticklabels=HIGHEST_DIVERGENCE,cmap=c,ax=ax,cbar=cbar)
    qm = ax.pcolormesh(T,cmap=c,edgecolors='none') 
    ax.set_frame_on(False)
    ax.set_yticks(list(range(len(HIGHEST_DIVERGENCE))),labels=HIGHEST_DIVERGENCE)

    return qm

fig, axs = plt.subplots(1,5,sharey=True,figsize=(15,12))
plot(DIR, ax=axs[0])
plot(DIR5, ax=axs[1])
plot(DIR4, ax=axs[2])
plot(DIR2, ax=axs[3])
qm = plot(DIR3, ax=axs[4])
for ax in axs:
    ax.set(xlabel="height in phylogeny")
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='x', rotation='default')
    ax.set_xticks([0,5,10])
    ax.set_xticklabels([0,5,10])
axs[0].set_title(r'RF_$\emptyset$ (but still DLC)')
axs[1].set_title(r'RE_$\emptyset$')
axs[2].set_title(r'IL_$\emptyset$')
axs[3].set_title('IL_ILC')
axs[4].set_title('RF_ILC')

fig.colorbar(qm, ax = axs[:]).set_label(label='number of base-pairs (normalized)', size=16)
fig.savefig('figures/rfam_matrix.pdf', bbox_inches='tight')
plt.show()
