import matplotlib.pylab as plt
from rnadist.utils import pairing_info, depthk_bps

DIR_HEUR_C2_IL = 'results/small_parsimony_results_random_input_median_heuristic/'
DIR_HEUR_C2_RF = 'results/small_parsimony_results_random_input_c2_rf_median_heuristic/'

N_random = 100


def num_loops(structure):


    def aux(i,j):
        if j <= i: 
            return 0

        b0 = depthk_bps(structure[i:j])
        b1 = depthk_bps(structure[i:j], depth=1)
        cnt = 0
        if len(b0)==1 and len(b1)>1:
            cnt += 1
        for k,l in b0:
            cnt += aux(k+1,l-1)

        return cnt

    return aux(0,len(structure))

x_values = []
y_values = []

for i in range(N_random):
    s_il = open('results/small_parsimony_results_random_input_median_heuristic/'+str(i)+'_fully_annotated_median_heuristic.tab').readlines()[0].rstrip('\n').split(':')[-1]
    s_rf = open('results/small_parsimony_results_random_input_c2_rf_median_heuristic/'+str(i)+'_fully_annotated_c2_rf_median_heuristic.tab').readlines()[0].rstrip('\n').split(':')[-1]
    print('predictions are the same', s_il==s_rf)
    x_values.append(s_il.count('('))
    y_values.append(s_rf.count('('))

plt.scatter(x_values, y_values)
plt.plot([4,7],[4,7],color='r')
plt.show()
