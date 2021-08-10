import matplotlib.pyplot as plt
import random
import numpy as np
import matplotlib as mpl


import matplotlib.font_manager
#print( matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf'))
#matplotlib.font_manager._rebuild()
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 24


dose_data = {}
dose_incorporated = {}

anneal_data = {}
anneal_incorporated = {}

pressures = ['3e-12','3e-11','3e-10','3e-9','3e-8']

times = [6,60,600,6000,60000,600000]

incorporated_mat = np.zeros((len(pressures),len(times)))

for i in range(len(pressures)):
    for j in range(len(times)):
        pressure = pressures[i]
        time = times[j]
        print (pressure)
        print (time)
        anneal_data = np.loadtxt('incorporation-A1e+12-200-runs-pressure%s-time%s.csv'%(pressure,time),delimiter = ',',skiprows=1)
        anneal_incorporated = {}

        for row in range(len(anneal_data[:,0])):
            anneal_incorporated[anneal_data[row,0]] = anneal_data[row,1:]

        # fuschle = { 3: np.array([0.29,0.71,0.0,0.0]),
        #             4: np.array([0.22,0.57,0.22,0.0]),
        #             5: np.array([0.15,0.55,0.30,0.0]),
        #             6: np.array([0.0,0.25,0.25,0.5])
        #             }

        proportion ={2: np.array([3.0/6.0,3.0/6.0]),
                          3: np.array([24.0/31.0,7.0/31.0])}

        combo_anneal = {}
        for value in proportion.keys():
            combo_anneal[value] = proportion[value][0]*anneal_incorporated[value] + proportion[value][1]*anneal_incorporated[value + 0.5]

        incorporated_mat[i,j] = combo_anneal[3][1]




print (incorporated_mat.transpose())



fig,ax = plt.subplots(figsize=(10,10))
# best so far is gnuplot with limits 0 to 1
plt.imshow(incorporated_mat.transpose(),origin='lower',cmap='Blues')
plt.clim(0,1)
cbar = plt.colorbar()
cbar.set_label('Proportion of P incorporation in 3 dimer window')
plt.plot([-0.5,4.5],[4.5,-0.5],c='k')
plt.text(2.0,2.0,'0.18 L',color = 'k')
plt.plot([-0.5,3.5],[3.5,-0.5],c='k')
plt.text(2.0,3.0,'1.8 L',color = 'k')
plt.plot([-0.5,2.5],[2.5,-0.5],c='k')
plt.text(2.0,4.0,'18 L',color = 'k')
plt.plot([-0.5,1.5],[1.5,-0.5],c='k')
plt.text(2.0,5.0,'180 L',color = 'k')
plt.plot([-0.5,0.5],[0.5,-0.5],c='k')
plt.text(2.0,1.0,'0.018 L',color = 'k')
plt.plot([-0.5,4.5],[5.5,0.5],c='k')
plt.text(2.0,0.0,'0.0018 L',color = 'k')
plt.plot([0.5,4.5],[5.5,1.5],c='k')
plt.plot([1.5,4.5],[5.5,2.5],c='k')
plt.plot([2.5,4.5],[5.5,3.5],c='k')
plt.plot([3.5,4.5],[5.5,4.5],c='k')

ax.set_xticks(range(len(pressures)))
ax.set_yticks(range(len(times)))
ax.set_xticklabels(pressures)
ax.set_yticklabels(times)

plt.xlabel('Dose Pressure (Torr)')
plt.ylabel('Dose Time (s)')


plt.tight_layout()
plt.savefig('heatmap.pdf')
plt.close()
#plt.show()
