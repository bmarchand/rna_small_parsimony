import seaborn as sns
import matplotlib.pylab as plt
import pandas as pd
import numpy as np
import json

with open('box_plot_values_random.json') as f:
    d = json.load(f)
with open('box_plot_values.json') as f:
    d2 = json.load(f)

data = {}
data['SP cost'] = []
data['method'] = []
data['distance'] = []
data2 = {}
data2['SP cost'] = []
data2['method'] = []
data2['distance'] = []

map_meth = { 0 : 'RF_NC', 1: 'IL_NC', 2:'RF_ILC',3:'IL_ILC',4:'RE',}
map_dist = { 0 : 'RF', 1: 'IL', 2:'RE'}

for k1 in d.keys():
    for k2, l in d[k1].items():
        for val in l:
            #data.append([float(val),map_meth[int(k1)],k2])
            data['SP cost'].append(float(val))
            data['method'].append(map_meth[int(k1)])
            data['distance'].append(map_dist[int(k2)])
for k1 in d2.keys():
    for k2, l in d2[k1].items():
        for val in l:
            #data.append([float(val),map_meth[int(k1)],k2])
            data2['SP cost'].append(float(val))
            data2['method'].append(map_meth[int(k1)])
            data2['distance'].append(map_dist[int(k2)])


df = pd.DataFrame(data, columns=['SP cost','method','distance'])
print(df.info())
df2 = pd.DataFrame(data2, columns=['SP cost','method','distance'])
print(df2.info())

fig, axs = plt.subplots(2)

sns.boxplot(data=df,ax=axs[0], hue="method",y=df["SP cost"],x='distance')
sns.boxplot(data=df2,ax=axs[1], hue="method",y=df2["SP cost"],x='distance')
axs[1].get_legend().remove()
fig.tight_layout()
plt.savefig('figures/box_plot_combined.pdf')
plt.show()
