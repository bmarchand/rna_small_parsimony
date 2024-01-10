from rnadist.utils import tree_tab_file_parse 
import os

numbps_values_per_height = {}
numbps_values_per_height_heuristic = {}

def height(node):
    if node.isleaf:
        return 0
    return min([height(child) for child in node.children])+1

DIR_FITCH = 'results/small_parsimony_results_random_input_fitch_rf//'
OUTNAME = 'figures/average_num_bps_height_random_structures.pdf'

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
        num = node.annotation.count('(')
        assert(node.annotation.count('(')==node.annotation.count(')'))
        try:
            numbps_values_per_height[h].append(num)
        except KeyError:
            numbps_values_per_height[h] = [num]

        for child in node.children:
            queue.append(child)


DIR_HEUR = 'results/small_parsimony_results_random_input_median_heuristic/'

# HEURISTIC
for filename in os.listdir(DIR_HEUR):
    print(filename)
    lines = open(DIR_HEUR+filename).readlines()

    phylo_T = tree_tab_file_parse(lines)

    print(phylo_T.label, phylo_T.annotation)

    queue = [phylo_T]

    while len(queue) > 0:
        node = queue.pop()
        h = height(node)
        num = node.annotation.count('(')
        assert(node.annotation.count('(')==node.annotation.count(')'))
        try:
            numbps_values_per_height_heuristic[h].append(num)
        except KeyError:
            numbps_values_per_height_heuristic[h] = [num]

        for child in node.children:
            queue.append(child)

import matplotlib.pylab as plt
import numpy as np

x_values = sorted(numbps_values_per_height.keys())

plt.violinplot([numbps_values_per_height[h] for h in x_values], positions=x_values,showmeans=False,showextrema=False)
medians = [np.median(numbps_values_per_height[h]) for h in x_values]
plt.plot(x_values, medians, 'o-', label='RF_fitch')
plt.boxplot([numbps_values_per_height[h] for h in x_values], positions=x_values,sym="")
plt.ylim([0,100])
plt.xlabel('height in phylogenetic tree')
plt.ylabel('average number of base pairs')
plt.savefig('figures/average_num_bps_height.pdf')
#plt.errorbar(x_values,[np.mean(numbps_values_per_height[h]) for h in x_values],yerr=[np.mean(numbps_values_per_height[h]) for h in x_values])

x_values_heuristic = sorted(numbps_values_per_height_heuristic.keys())
plt.violinplot([numbps_values_per_height_heuristic[h] for h in x_values_heuristic], positions=x_values_heuristic,showmeans=False,showextrema=False)
medians_heuristic = [np.median(numbps_values_per_height_heuristic[h]) for h in x_values_heuristic]
plt.plot(x_values_heuristic, medians_heuristic, 'o-', label='c2_il_median_heuristic')
plt.boxplot([numbps_values_per_height_heuristic[h] for h in x_values_heuristic], positions=x_values_heuristic,sym="")
plt.ylim([0,100])
plt.xlabel('height in phylogenetic tree')
plt.ylabel('average number of base pairs')
plt.legend()
plt.savefig(OUTNAME)
