import os
from rnadist.utils import tree_tab_file_parse
import json

FILTERED_RFAM = [line.rstrip('\n') for line in open('FILTERED_RFAM').readlines()]
#VERY_FILTERED_RFAM = [line.rstrip('\n') for line in open('VERY_FILTERED_RFAM').readlines()]
import json
with open('divergence.json') as f:
    divergence = json.load(f)

FILTERED_RFAM = list(sorted(FILTERED_RFAM, key=lambda x: -divergence[x]))[:20]
print(len(FILTERED_RFAM),'families')


DIR1 = 'results/small_parsimony_results_rfam_fitch_rf/'
DIR2 = 'results/small_parsimony_results_rfam_unconstrained_il_median_heuristic/'
DIR3 = 'results/small_parsimony_results_rfam_c2_rf_median_heuristic/'
DIR4 = 'results/small_parsimony_results_rfam_median_heuristic/'
DIR5 = 'results/small_parsimony_results_rfam_re_distance_heuristic/'

from rnadist.utils import RFdistance, IL_distance
from rnadist.RE_distance import RE_distance

with open('rf_distance.json') as f:
    c_rf = json.load(f)
with open('il_distance.json') as f:
    c_il = json.load(f)
with open('re_distance.json') as f:
    c_re = json.load(f)

def compute_RF(s1,s2):
    if s1+'-'+s2 in c_rf.keys():
        return c_rf[s1+'-'+s2]
    
    c_rf[s1+'-'+s2] = RFdistance(s1,s2)
    return c_rf[s1+'-'+s2]

def compute_IL(s1,s2):
    if s1+'-'+s2 in c_il.keys():
        return c_il[s1+'-'+s2]
    
    c_il[s1+'-'+s2] = IL_distance(s1,s2)
    return c_il[s1+'-'+s2]

def compute_RE(s1,s2):
    if s1+'-'+s2 in c_re.keys():
        return c_re[s1+'-'+s2]
    
    c_re[s1+'-'+s2] = RE_distance(s1,s2)
    return c_re[s1+'-'+s2]


def values(DIR, distance):
    vals = []
    for fname in os.listdir(DIR):
        family = fname.split('_')[0]
        if family not in FILTERED_RFAM:
            continue
        phylo_T = tree_tab_file_parse(open(DIR+fname).readlines())

        queue = [(phylo_T,c) for c in phylo_T.children]

        vals_fname = []
        while len(queue) > 0:
            u,v = queue.pop()
            if u.annotation==v.annotation:
                dist = 0
            else:
                dist =  distance(u.annotation,v.annotation)
            vals_fname.append(dist)

            for w in v.children:
                queue.append((v,w))
        
        vals.append(np.mean(vals_fname))

    return vals

import numpy as np
data = np.zeros((5,3))

values_dict = {}

for i, DIR in enumerate([DIR1,DIR2,DIR3,DIR4,DIR5]):
    values_dict[i] = {}
    for j, distance in enumerate([compute_RF,compute_IL,compute_RE]):
        vals = values(DIR, distance)
        print(vals)
        values_dict[i][j] = vals
#        data[i,j] = np.mean(values(DIR, distance))

import json
with open('box_plot_values.json','w') as f:
    json.dump(values_dict, f)

with open('rf_distance.json','w') as f:
    json.dump(c_rf, f)
with open('il_distance.json','w') as f:
    json.dump(c_il, f)
with open('re_distance.json','w') as f:
    json.dump(c_re, f)


#import matplotlib.pylab as plt
#import seaborn as sns
#sns.heatmap(data, yticklabels=['RF_ILC','RE','IL_ILC','RF_NC','IL_NC'],xticklabels=['RF','IL','RE'])
#plt.show()
