import os

FAMILIES = [fam.split('/')[-1].split('.')[0] for fam in os.listdir('resources/tree_files/')]

# Filtering families to retain only interesting ones, i.e:
#       1. For which there is agreement between the seed file and tree file compositions
#       2. With non empty structures at the leaves
not_nice = []
for family in FAMILIES:
    max_num_bps = 0
    try:
        if len(open('results/small_parsimony_input_rfam/'+family+'_leaf_annotated.tab').readlines()) > 10000000000:
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

divergence = {}

from rnadist.utils import IL_distance

for family in NICE_FAMILIES:
    print('computing divergence of',family)
    structures = []
    for line in open('results/small_parsimony_input_rfam/'+family+'_leaf_annotated.tab').readlines():
        if line.find(':') >= 0:
            struct = line.split(' ')[-1].rstrip('\n')
            structures.append(struct)

    sc = 0
    for i, s1 in enumerate(structures):
        for j, s2 in enumerate(structures):
            if i<j:
                sc += IL_distance(s1,s2)

    divergence[family] = float(sc) / float(len(structures)*(len(structures)-1)*0.5)
    print(divergence[family])

for k,v in sorted(divergence.items(),key=lambda x:x[1]):
    print(k,v)

import json
with open('divergence.json','w') as f:
    json.dump(divergence, f)
