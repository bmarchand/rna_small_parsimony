import matplotlib.pylab as plt
import os

FAMILIES = [fam.split('/')[-1].split('.')[0] for fam in os.listdir('resources/tree_files/')]

# Filtering families to retain only interesting ones, i.e:
#       1. For which there is agreement between the seed file and tree file compositions
#       2. With non empty structures at the leaves
not_nice = []
for family in FAMILIES:
    max_num_bps = 0
    try:
        if len(open('results/small_parsimony_input_rfam/'+family+'_leaf_annotated.tab').readlines()) > 100:
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

DIR1 = 'results/small_parsimony_results_rfam_fitch_rf/'
DIR2 = 'results/small_parsimony_results_rfam_median_heuristic/'
DIR3 = 'results/small_parsimony_results_rfam_c2_rf_median_heuristic/'
DIR4 = 'results/small_parsimony_results_rfam_unconstrained_il_median_heuristic/'

def value_dict(DIR):
    d = {}
    for fname in os.listdir(DIR):
        family = fname.split('_')[0]

        root_struct = open(DIR+fname).readlines()[0].split(' ')[-1].rstrip('\n')

        d[family] = root_struct.count('(')

    return d

d1 = value_dict(DIR1)
d2 = value_dict(DIR2)
d3 = value_dict(DIR3)
#d4 = value_dict(DIR4)

fig, axs = plt.subplots(3,3,sharex=True,sharey=True)

for i in range(1,3,1):
    for j in range(i):
        axs[i,j].axis('off')

axs[0,0].scatter([d2[f] for f in NICE_FAMILIES],[d1[f] for f in NICE_FAMILIES],s=6, c='black')
axs[0,0].set_xlabel('IL with IIL constraint')
axs[0,0].set_ylabel('unconstrained RF')
axs[0,0].plot([0,50],[0,50],color='r')
axs[0,0].xaxis.set_label_position('top')

axs[0,1].scatter([d3[f] for f in NICE_FAMILIES],[d1[f] for f in NICE_FAMILIES],s=6, c='black')
axs[0,1].plot([0,50],[0,50],color='r')
axs[0,1].set_xlabel('RF with IIL constraint')
axs[0,1].xaxis.set_label_position('top')

axs[1,1].scatter([d3[f] for f in NICE_FAMILIES],[d2[f] for f in NICE_FAMILIES],s=6, c='black')
axs[1,1].plot([0,50],[0,50],color='r')
axs[1,1].set_ylabel('IL with IIL constraint')

fig.suptitle('Comparison of the number of base-pairs\n predicted at the root of the phylogeny for each RFAM family')
plt.show()
