import os
import json

with open('id2family.json') as f:
    id2family = json.load(f)

family_per_clan = {}
for line in open('resources/Rfam.clanin').readlines():
    clan = line.split('\t')[0]
    family_per_clan[clan] = []
    for ident in line.rstrip('\n').split('\t')[1:]:
        family_per_clan[clan].append(ident)

for clan, idents in family_per_clan.items():
    print(clan)
    for ident in idents:
        family = id2family[ident]

        os.system('cp ../rf_small_parsimony/results/small_parsimony_results_rfam_median_heuristic/'+family+'_fully_annotated_median_heuristic.tab resources/annotated_family_trees/'+clan+'/')
        os.system('cp ../rf_small_parsimony/results/small_parsimony_results_rfam_c2_rf_median_heuristic/'+family+'_fully_annotated_c2_rf_median_heuristic.tab resources/annotated_family_trees/'+clan+'/')
        os.system('cp ../rf_small_parsimony/results/small_parsimony_results_rfam_fitch_rf/'+family+'_fully_annotated_fitch_rf.tab resources/annotated_family_trees/'+clan+'/')
