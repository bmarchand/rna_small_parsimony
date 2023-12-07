from rnadist.utils import tree_tab_file_parse 
import os

degree_values_per_height = {}

def height(node):
    if node.isleaf:
        return 0
    return min([height(child) for child in node.children])+1

def multiloop_degrees(structure):



for filename in os.listdir('results/small_parsimony_results/'):
    print(filename)
    lines = open('results/small_parsimony_results/'+filename).readlines()

    phylo_T = tree_tab_file_parse(lines)

    print(phylo_T.label, phylo_T.annotation)
