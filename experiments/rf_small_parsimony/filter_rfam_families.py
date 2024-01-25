import os
import json

FAMILIES = [fam.split('/')[-1].split('.')[0] for fam in os.listdir('resources/tree_files/')]

with open('divergence.json') as f:
    divergence = json.load(f)

# Filtering families to retain only interesting ones, i.e:
#       1. For which there is agreement between the seed file and tree file compositions
#       2. With non empty structures at the leaves
phylo_too_large = []
rna_way_too_large = []
rna_too_large = []
no_structure = []
not_found = []
for family in FAMILIES:
    max_num_bps = 0
    try:
        if len(open('results/small_parsimony_input_rfam/'+family+'_leaf_annotated.tab').readlines()) > 1000000000000000:
            phylo_too_large.append(family)
        for line in open('results/small_parsimony_input_rfam/'+family+'_leaf_annotated.tab').readlines():
            if line.find(':') >= 0:
                struct = line.split(' ')[-1].rstrip('\n')
                max_num_bps = max(max_num_bps, struct.count('('))
        if max_num_bps==0:
            no_structure.append(family)
        if len(struct) > 30:
            rna_too_large.append(family)
        if len(struct) > 100:
            rna_way_too_large.append(family)
    except FileNotFoundError:
        not_found.append(family)

FILTERED_RFAM = [fam for fam in FAMILIES if fam not in rna_way_too_large and fam not in no_structure and fam not in not_found]
FILTERED_RFAM = [fam for fam in FILTERED_RFAM if divergence[fam] > 0]
VERY_FILTERED_RFAM = [fam for fam in FAMILIES if fam not in rna_too_large and fam not in phylo_too_large and fam not in no_structure and fam not in not_found]
VERY_FILTERED_RFAM = [fam for fam in VERY_FILTERED_RFAM if divergence[fam] > 0]

print('FILTERED_RFAM:', len(FILTERED_RFAM),'families')
print('VERY_FILTERED_RFAM:', len(VERY_FILTERED_RFAM),'families')

with open('FILTERED_RFAM','w') as f:
    for fam in FILTERED_RFAM:
        print(fam, file=f)

with open('VERY_FILTERED_RFAM','w') as f:
    for fam in VERY_FILTERED_RFAM:
        print(fam, file=f)
