import os

FILTERED_RFAM = [line.rstrip('\n') for line in open('FILTERED_RFAM').readlines()]
VERY_FILTERED_RFAM = [line.rstrip('\n') for line in open('VERY_FILTERED_RFAM').readlines()]

divergence = {}

from rnadist.utils import IL_distance

for family in FILTERED_RFAM:
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
