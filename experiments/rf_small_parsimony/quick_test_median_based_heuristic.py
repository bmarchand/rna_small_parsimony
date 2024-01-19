from rnadist.utils import tree_tab_file_parse, median_based_heuristic

# parsing input tree file
phylo_T = tree_tab_file_parse(open('results/small_parsimony_input_rfam/RF00004_leaf_annotated.tab').readlines())

def pprint(node, depth):
    print('\t'*depth, node.label)
    for child in node.children:
        pprint(child, depth+1)

pprint(phylo_T,0)

# making dict: leaf_label -> str
str_dict = {}
for line in open('results/small_parsimony_input_rfam/RF00004_leaf_annotated.tab').readlines():
    if line.find(':') >= 0:
        acc = line.split(' ')[-3].lstrip('\t')
        struct = line.split(' ')[-1].rstrip('\n')
        str_dict[acc] = struct

str_dict = median_based_heuristic(phylo_T, str_dict)
