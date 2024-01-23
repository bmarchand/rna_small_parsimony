from rnadist.utils import tree_tab_file_parse, pairing_info, depthk_bps 
import sys
import os

FAMILIES = [fam.split('/')[-1].split('.')[0] for fam in os.listdir('resources/tree_files/')]

# Filtering families to retain only interesting ones, i.e:
#       1. For which there is agreement between the seed file and tree file compositions
#       2. With non empty structures at the leaves
not_nice = []
for family in FAMILIES:
    max_num_bps = 0
    try:
        if len(open('results/small_parsimony_input_rfam/'+family+'_leaf_annotated.tab').readlines()) > 100:
            not_nice.append(family)
            continue
        for line in open('results/small_parsimony_input_rfam/'+family+'_leaf_annotated.tab').readlines():
#            if line.find('UNKNOWN') >= 0:           
#                not_nice.append(family)
#                break       
            if line.find(':') >= 0:
                struct = line.split(' ')[-1].rstrip('\n')
                max_num_bps = max(max_num_bps, struct.count('('))
        if max_num_bps==0 or len(struct) > 100:
            not_nice.append(family)
    except FileNotFoundError:
        not_nice.append(family)

NICE_FAMILIES = [fam for fam in FAMILIES if fam not in not_nice]


if len(sys.argv)==1: 
    dataset = 'RFAM'
    #dataset = 'random'
    height_method = 'max'
    #height_method = 'min'
    
    metric = 'num_bps'
    #metric = 'num_loops'
else:
    height_method = sys.argv[-1]
    metric = sys.argv[-2]
    dataset = sys.argv[-3]

if dataset=='RFAM':
    DIR_FITCH = 'results/small_parsimony_results_rfam_fitch_rf/'
    DIR_HEUR_C2_IL = 'results/small_parsimony_results_rfam_median_heuristic/'
    DIR_HEUR_C2_RF = 'results/small_parsimony_results_rfam_c2_rf_median_heuristic/'
    DIR_HEUR_UNC_IL = 'results/small_parsimony_results_rfam_unconstrained_il_median_heuristic/'

    if height_method=='max':
        OUTNAME = 'figures/average_'+metric+'_maxheight_rfam.pdf'
    if height_method=='min':
        OUTNAME = 'figures/average_'+metric+'_height_rfam.pdf'

if dataset=='random':
    DIR_FITCH = 'results/small_parsimony_results_random_input_fitch_rf/'
    DIR_HEUR_C2_IL = 'results/small_parsimony_results_random_input_median_heuristic/'
    DIR_HEUR_C2_RF = 'results/small_parsimony_results_random_input_c2_rf_median_heuristic/'
    DIR_HEUR_UNC_IL = 'results/small_parsimony_results_random_input_unconstrained_il_median_heuristic/'

    if height_method=='max':
        OUTNAME = 'figures/average_'+metric+'_maxheight_random_structures.pdf'
    if height_method=='min':
        OUTNAME = 'figures/average_'+metric+'_height_random_structures.pdf'

def num_loops(structure):

    #d, _ = pairing_info(structure

    def aux(i,j, depth=0):
        if j <= i+1: 
            return 0

        print('\t'*depth+structure[i:j])
        l0 = depthk_bps(structure[i:j])
        print('\t'*depth+'l0',l0)
    
        cnt = 0
        if len(l0) > 1:
            cnt += 1
        for k,l in l0:
            cnt += aux(i+k+1,i+l-1, depth=depth+1)

        return cnt

    return aux(0,len(structure))


def height(node):
    if node.isleaf:
        return 0
    if height_method=='max':
        return max([height(child) for child in node.children])+1
    if height_method=='min':
        return min([height(child) for child in node.children])+1

def values_per_height_dict(DIR):
    d = {}

    for filename in os.listdir(DIR):
        print(filename)
        if  dataset=='RFAM' and filename.split('_')[0] not in NICE_FAMILIES:
            continue
        lines = open(DIR+filename).readlines()
    
        phylo_T = tree_tab_file_parse(lines)
    
        queue = [phylo_T]

        d_filename = {}
    
        while len(queue) > 0:
            node = queue.pop()
            h = height(node)
            if metric=='num_bps':
                num = node.annotation.count('(')
            if metric=='num_loops':
                num = num_loops(node.annotation)
            assert(node.annotation.count('(')==node.annotation.count(')'))
            try:
                d_filename[h].append(num)
            except KeyError:
                d_filename[h] = [num]
    
            for child in node.children:
                queue.append(child)

        for h, val in d_filename.items():
            try:
                d[h].append(max(val))
            except KeyError:
                d[h] = [max(val)]

    return d

# FITCH
numbps_values_per_height = values_per_height_dict(DIR_FITCH)

# HEURISTIC C2_IL
numbps_values_per_height_c2_il_heuristic = values_per_height_dict(DIR_HEUR_C2_IL)

# HEURISTIC C2_RF
numbps_values_per_height_c2_rf_heuristic = values_per_height_dict(DIR_HEUR_C2_RF)

# HEURISTIC UNC_IL
numbps_values_per_height_unconstrained_il_heuristic = values_per_height_dict(DIR_HEUR_UNC_IL)

import matplotlib.pylab as plt
import seaborn as sns
import numpy as np

sns.set_palette('colorblind')

def plot(d, label):
    x_values = sorted(d.keys())
    medians = [np.mean(d[h]) for h in x_values]
    stds = [np.std(d[h])/np.sqrt(len(d[h])) for h in x_values]
    plt.errorbar(x_values, 
                 medians, 
                 yerr=stds, 
                 fmt='o-',
                 label=label,
                 markersize=3,
                 capsize=2)

plot(numbps_values_per_height, 'fitch_RF')
plot(numbps_values_per_height_c2_il_heuristic, 'c2_il_heuristic')
plot(numbps_values_per_height_c2_rf_heuristic, 'c2_rf_heuristic')
plot(numbps_values_per_height_unconstrained_il_heuristic, 'unc_heuristic')

fontsize = 15

plt.xlabel('height in phylogenetic tree',fontsize=fontsize)
#plt.title(dataset+' x '+metric+' x '+height_method)


if metric=='num_bps':
    plt.ylabel('average number of base pairs',fontsize=fontsize)
if metric=='num_loops':
    plt.ylabel('average number of multi-loops',fontsize=fontsize)


plt.title('Average maximum number of base-pairs\n in predicted ancestral structures (RANDOM dataset)')
plt.legend()
plt.savefig(OUTNAME)
#plt.show()
