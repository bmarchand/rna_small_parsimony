import os
import numpy as np


FILTERED_RFAM = [line.rstrip('\n') for line in open('FILTERED_RFAM').readlines()]
VERY_FILTERED_RFAM = [line.rstrip('\n') for line in open('VERY_FILTERED_RFAM').readlines()]

values = {'rf_nc':[],'il_nc':[],'il_ilc':[],'rf_ilc':[],'re':[]}

for fname in os.listdir('benchmarks/'):
    if fname.find('rfam') >= 0:

        family = fname.split('_')[-2]
        if family not in FILTERED_RFAM:
            continue

        if fname.startswith('c2_rf_median'):
            values['rf_ilc'].append(float(open('benchmarks/'+fname).readlines()[1].split('\t')[0]))
        if fname.startswith('median'):
            values['il_ilc'].append(float(open('benchmarks/'+fname).readlines()[1].split('\t')[0]))
        if fname.startswith('unconstrained'):
            values['il_nc'].append(float(open('benchmarks/'+fname).readlines()[1].split('\t')[0]))
        if fname.startswith('fitch'):
            values['rf_nc'].append(float(open('benchmarks/'+fname).readlines()[1].split('\t')[0]))
        if fname.startswith('re'):
            values['re'].append(float(open('benchmarks/'+fname).readlines()[1].split('\t')[0]))
        
for k, v in values.items():
    print(k, np.mean(v))
