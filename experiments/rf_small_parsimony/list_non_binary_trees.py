from rnadist.utils import tree_tab_file_parse
import os

binary = 0
non_binary = 0
some_unary = set([])
max_degrees = []

for filename in os.listdir('results/unannotated_tree_files/'):
    tree_lines = open('results/unannotated_tree_files/'+filename).readlines()
    phylo_T = tree_tab_file_parse(tree_lines)

    max_degree = 0

    queue = [phylo_T]

    while len(queue) > 0:
        node = queue.pop()

        max_degree = max(max_degree, len(node.children))
        if len(node.children)==1:
            print('degree 1', filename, node.label)
            some_unary.add(filename)


        for child in node.children:
            queue.append(child)

    if max_degree <= 2:
        binary += 1
        # print(filename, 'max degree', max_degree)
    else:
        max_degrees.append(max_degree)
        non_binary += 1

import numpy as np
print('number of binary trees:', binary)
print('number of non binary trees:', non_binary,'avg max degree for them:', np.mean(max_degrees))
print('number of trees containing some unary nodes:', len(some_unary))
