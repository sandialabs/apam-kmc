import matplotlib.pyplot as plt
import random
import numpy as np
import matplotlib as mpl


import matplotlib.font_manager
#print( matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext='ttf'))
#matplotlib.font_manager._rebuild()
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 28


dose_data = {}
dose_incorporated = {}

anneal_data = {}
anneal_incorporated = {}

times = [9,90,900,9000,90000]
pressures = ['4e-11','4e-10','4e-9','4e-8']

incorporated_mat = np.zeros((len(pressures),len(times)))

for i in range(len(pressures)):
    for j in range(len(times)):
        pressure = pressures[i]
        time = times[j]
        anneal_data = np.loadtxt('incorporation-200-runs-single-acceptor-P-%s-t-%s.csv'%(pressure,time),delimiter = ',',skiprows=1)
        anneal_incorporated = {}
        anneal_incorporated[3] = anneal_data[1:]

        # # fuschle = { 3: np.array([0.29,0.71,0.0,0.0]),
        # #             4: np.array([0.22,0.57,0.22,0.0]),
        # #             5: np.array([0.15,0.55,0.30,0.0]),
        # #             6: np.array([0.0,0.25,0.25,0.5])
        # #             }
        #
        # proportion ={2: np.array([3.0/6.0,3.0/6.0]),
        #                   3: np.array([24.0/31.0,7.0/31.0])}
        #
        # combo_anneal = {}
        # for value in proportion.keys():
        #     combo_anneal[value] = proportion[value][0]*anneal_incorporated[value] + proportion[value][1]*anneal_incorporated[value + 0.5]

        incorporated_mat[i,j] = anneal_incorporated[3][1]




print (incorporated_mat.transpose())


# cmap_schemes = ['viridis', 'plasma', 'inferno', 'magma', 'cividis','Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
#             'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
#             'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
#             'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
#             'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
#             'hot', 'afmhot', 'gist_heat', 'copper',
#             'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
#             'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
#             'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']

cmap_schemes = [ 'Blues']


for cmap in cmap_schemes:
    fig,ax = plt.subplots(figsize=(10,10))
    # best so far is gnuplot with limits 0 to 1
    plt.imshow(incorporated_mat.transpose(),origin='lower',cmap=cmap)
    plt.clim(0,1)
    cbar = plt.colorbar()
    cbar.set_label('Probability of Al incorporation')
    plt.plot([-0.5,3.5],[4.5,0.5],c='k')
    plt.text(1.0,2.0,'0.36 L',color = 'w')
    plt.plot([-0.5,3.5],[3.5,-0.5],c='k')
    plt.text(1.0,3.0,'3.6 L',color = 'w')
    plt.plot([-0.5,2.5],[2.5,-0.5],c='k')
    plt.text(1.0,4.0,'36 L',color = 'w')
    plt.plot([-0.5,1.5],[1.5,-0.5],c='k')
    plt.plot([-0.5,0.5],[0.5,-0.5],c='k')
    plt.text(1.0,1.0,'0.036 L',color = 'k')
    plt.plot([0.5,3.5],[4.5,1.5],c='k')
    plt.text(1.0,0.0,'0.0036 L',color = 'k')
    plt.plot([1.5,3.5],[4.5,2.5],c='k')
    plt.plot([2.5,3.5],[4.5,3.5],c='k')

    ax.set_xticks(range(len(pressures)))
    ax.set_yticks(range(len(times)))
    ax.set_xticklabels(pressures)
    ax.set_yticklabels(times)

    plt.xlabel('Dose Pressure (Torr)')
    plt.ylabel('Dose Time (s)')


    plt.tight_layout()
    plt.savefig('./heatmaps/%s_heatmap.pdf'%cmap)
    #plt.close()
    plt.show()
