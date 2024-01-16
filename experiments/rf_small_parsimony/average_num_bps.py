from rnadist.utils import tree_tab_file_parse, pairing_info 
import os

numbps_values_per_height = {}
numbps_values_per_height_c2_il_heuristic = {}
numbps_values_per_height_c2_rf_heuristic = {}

#dataset = 'RFAM'
dataset = 'random'
#height_method = 'max'
height_method = 'min'

#metric = 'num_bps'
metric = 'num_loops'

if dataset=='RFAM':
    DIR_FITCH = 'results/small_parsimony_results_rfam_fitch_rf/'
    DIR_HEUR_C2_IL = 'results/small_parsimony_results_rfam_median_heuristic/'
    DIR_HEUR_C2_RF = 'results/small_parsimony_results_rfam_c2_rf_median_heuristic/'
    if height_method=='max':
        OUTNAME = 'figures/average_'+metric+'_maxheight_rfam.pdf'
    if height_method=='min':
        OUTNAME = 'figures/average_'+metric+'_height_rfam.pdf'

if dataset=='random':
    DIR_FITCH = 'results/small_parsimony_results_random_input_fitch_rf/'
    DIR_HEUR_C2_IL = 'results/small_parsimony_results_random_input_median_heuristic/'
    DIR_HEUR_C2_RF = 'results/small_parsimony_results_random_input_c2_rf_median_heuristic/'
    if height_method=='max':
        OUTNAME = 'figures/average_'+metric+'_maxheight_random_structures.pdf'
    if height_method=='min':
        OUTNAME = 'figures/average_'+metric+'_height_random_structures.pdf'

def num_loops(structure):

    d, _ = pairing_info(structure)

    def aux(i,j):
        if j <= i: 
            return 0
        if i in d.keys():
            return 1+aux(i+1,d[i]-1)+aux(d[i]+1,j)
        else:
            return aux(i+1,j)

    return aux(0,len(structure)-1)


def height(node):
    if node.isleaf:
        return 0
    if height_method=='max':
        return max([height(child) for child in node.children])+1
    if height_method=='min':
        return min([height(child) for child in node.children])+1


# FITCH
for filename in os.listdir(DIR_FITCH):
    print(filename)
    lines = open(DIR_FITCH+filename).readlines()

    phylo_T = tree_tab_file_parse(lines)

    print(phylo_T.label, phylo_T.annotation)

    queue = [phylo_T]

    while len(queue) > 0:
        node = queue.pop()
        h = height(node)
        if metric=='num_bps':
            num = node.annotation.count('(')
        if metric=='num_loops':
            num = num_loops(node.annotation)
        assert(node.annotation.count('(')==node.annotation.count(')'))
        try:
            numbps_values_per_height[h].append(num)
        except KeyError:
            numbps_values_per_height[h] = [num]

        for child in node.children:
            queue.append(child)


# HEURISTIC C2_IL
for filename in os.listdir(DIR_HEUR_C2_IL):
    print(filename)
    lines = open(DIR_HEUR_C2_IL+filename).readlines()

    phylo_T = tree_tab_file_parse(lines)

    print(phylo_T.label, phylo_T.annotation)

    queue = [phylo_T]

    while len(queue) > 0:
        node = queue.pop()
        h = height(node)
        num = node.annotation.count('(')
        assert(node.annotation.count('(')==node.annotation.count(')'))
        try:
            numbps_values_per_height_c2_il_heuristic[h].append(num)
        except KeyError:
            numbps_values_per_height_c2_il_heuristic[h] = [num]

        for child in node.children:
            queue.append(child)

# HEURISTIC C2_RF
for filename in os.listdir(DIR_HEUR_C2_RF):
    print(filename)
    lines = open(DIR_HEUR_C2_RF+filename).readlines()

    phylo_T = tree_tab_file_parse(lines)

    print(phylo_T.label, phylo_T.annotation)

    queue = [phylo_T]

    while len(queue) > 0:
        node = queue.pop()
        h = height(node)
        num = node.annotation.count('(')
        assert(node.annotation.count('(')==node.annotation.count(')'))
        try:
            numbps_values_per_height_c2_rf_heuristic[h].append(num)
        except KeyError:
            numbps_values_per_height_c2_rf_heuristic[h] = [num]

        for child in node.children:
            queue.append(child)


import matplotlib.pylab as plt
import numpy as np

x_values = sorted(numbps_values_per_height.keys())
#plt.violinplot([numbps_values_per_height[h] for h in x_values], positions=x_values,showmeans=False,showextrema=False)
medians = [np.median(numbps_values_per_height[h]) for h in x_values]
plt.plot(x_values, medians, 'o-', label='RF_fitch')
plt.boxplot([numbps_values_per_height[h] for h in x_values], positions=x_values,sym="")

x_values_heuristic = sorted(numbps_values_per_height_c2_il_heuristic.keys())
#plt.violinplot([numbps_values_per_height_c2_il_heuristic[h] for h in x_values_heuristic], positions=x_values_heuristic,showmeans=False,showextrema=False)
medians_heuristic = [np.median(numbps_values_per_height_c2_il_heuristic[h]) for h in x_values_heuristic]
plt.plot(x_values_heuristic, medians_heuristic, 'o-', label='c2_il_heuristic')
plt.boxplot([numbps_values_per_height_c2_il_heuristic[h] for h in x_values_heuristic], positions=x_values_heuristic,sym="")

x_values_heuristic = sorted(numbps_values_per_height_c2_rf_heuristic.keys())
#plt.violinplot([numbps_values_per_height_c2_il_heuristic[h] for h in x_values_heuristic], positions=x_values_heuristic,showmeans=False,showextrema=False)
medians_heuristic = [np.median(numbps_values_per_height_c2_rf_heuristic[h]) for h in x_values_heuristic]
plt.plot(x_values_heuristic, medians_heuristic, 'o-', label='c2_rf_heuristic')
plt.boxplot([numbps_values_per_height_c2_il_heuristic[h] for h in x_values_heuristic], positions=x_values_heuristic,sym="")


plt.ylim([0,100])
plt.xlabel('height in phylogenetic tree')


if metric=='num_bps':
    plt.ylabel('average number of base pairs')
if metric=='num_loops':
    plt.ylabel('average number of loops')


plt.legend()
plt.savefig(OUTNAME)
plt.show()
