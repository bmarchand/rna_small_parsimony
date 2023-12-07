from rnadist.utils import tree_tab_file_parse 
import os

numbps_values_per_height = {}

def height(node):
    if node.isleaf:
        return 0
    return min([height(child) for child in node.children])+1

for filename in os.listdir('results/small_parsimony_results/'):
    print(filename)
    lines = open('results/small_parsimony_results/'+filename).readlines()

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

import matplotlib.pylab as plt
import numpy as np

x_values = sorted(numbps_values_per_height.keys())

print(numbps_values_per_height)

plt.violinplot([numbps_values_per_height[h] for h in x_values], positions=x_values,showmeans=False,showextrema=False)
medians = [np.median(numbps_values_per_height[h]) for h in x_values]
plt.plot(x_values, medians, 'o-')
plt.boxplot([numbps_values_per_height[h] for h in x_values], positions=x_values,sym="")
plt.ylim([0,100])
plt.xlabel('height in phylogenetic tree')
plt.ylabel('average number of base pairs')
plt.savefig('figures/average_num_bps_height.pdf')
#plt.errorbar(x_values,[np.mean(numbps_values_per_height[h]) for h in x_values],yerr=[np.mean(numbps_values_per_height[h]) for h in x_values])
