import seaborn as sns
import matplotlib.pylab as plt
import pandas as pd
import numpy as np
import json

with open('box_plot_values.json') as f:
    d = json.load(f)

data = {}
data['SP score'] = []
data['method'] = []
data['distance'] = []

map_meth = { 0 : 'RF_ILC', 1: 'RE', 2:'IL_ILC',3:'RF_NC',4:'IL_NC',}
map_dist = { 0 : 'RF', 1: 'IL', 2:'RE'}

for k1 in d.keys():
    for k2, l in d[k1].items():
        for val in l:
            #data.append([float(val),map_meth[int(k1)],k2])
            data['SP score'].append(float(val))
            data['method'].append(map_meth[int(k1)])
            data['distance'].append(map_dist[int(k2)])

np_data = np.array(data)
print(np_data.shape)

df = pd.DataFrame(data, columns=['SP score','method','distance'])
print(df.info())

sns.boxplot(data=df, hue="method",y=df["SP score"],x='distance')
plt.savefig('figures/box_plot.pdf')
plt.show()
