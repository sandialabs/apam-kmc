#Quinn Campbell qcampbe@sandia.gov

from KMCLib import *
import matplotlib.pyplot as plt
import random
import numpy as np
import matplotlib as mpl
import time
import math
import json

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Avenir'

dimer_height = 0.766
dimer_width = 0.3874

def num_avail_dimers(width,height):
    #for a given width and height, determines how many dimer sites
    #are theoretically available

    num_dimers_wide = math.ceil(width/dimer_width)
    num_dimers_tall = math.ceil(height/dimer_height)
    return num_dimers_wide*num_dimers_tall*2.0

def num_dimers_xy(width,height):
    #for a given width and height, determines how many dimers
    #are generated
    num_dimers_wide = math.ceil(width/dimer_width)
    num_dimers_tall = math.ceil(height/dimer_height)
    return num_dimers_wide+4, num_dimers_tall+3

nruns = 200

inserted= {}
inserted_loc ={}
wall_time = {}

distances = []
# here asssuming a patch of size patch_size nm x patch_size nm
# the config file will take care of converting that to a window
patch_size = [10.0]
for dimer_window in patch_size:
    inserted[dimer_window] = []
    start_time = time.time()
    time_to_check_after = 600.0  #In units of s
    inserted_loc[dimer_window] = np.zeros(num_dimers_xy(dimer_window,dimer_window))
    for i in range(nruns):
        time_to_check_after = 600.0  #In units of s
        num_incorporated = 0
        incorp_sites = []
        interactions  = KMCInteractionsFromScript('diborane-rxn.py')
        control_parameters = KMCControlParameters(number_of_steps=int(8*dimer_window*dimer_window),
                                                    dump_interval=int(8),seed=random.randint(1,2e9))
        configuration = KMCConfigurationFromScript('%s-%sconfig.py'%(dimer_window,dimer_window))
        model = KMCLatticeModel(configuration,interactions)
        try:
            model.run(control_parameters,trajectory_filename='%s-%spatch.py'%(dimer_window,dimer_window))
        except:
            print ( ' No more options available')


        # trajfile='%s-%spatch.py'%(dimer_window,dimer_window)
        # global_dict = {}
        # local_dict  = {}
        # execfile(trajfile, global_dict, local_dict)
        # times = local_dict['times']
        # sites = local_dict["sites"]
        # elem = local_dict["types"]
        interactions = KMCInteractionsFromScript('diborane-rxn-highT.py')
        control_parameters = KMCControlParameters(number_of_steps=int(40.0*dimer_window*dimer_window),
                                                    dump_interval=int(40.0),seed=random.randint(1,2e9))
        configuration = KMCConfigurationFromScript('%s-%sconfig_highT.py'%(dimer_window,dimer_window))
        model = KMCLatticeModel(configuration,interactions)
        try:
            model.run(control_parameters,trajectory_filename='%s-%spatch-T.py'%(dimer_window,dimer_window))
        except:
            print ( ' No more options available')

        trajfile='%s-%spatch-T.py'%(dimer_window,dimer_window)
        global_dict = {}
        local_dict  = {}
        execfile(trajfile, global_dict, local_dict)
        times = local_dict['times']
        sites = local_dict["sites"]
        elem = local_dict["types"]

        time_to_check_after = 60.0  #In units of s



        final_index = -1
        # for j in range(len(times)):
        #     print times[j]
        #     if times[j] > time_to_check_after:
        #         final_index = j-1
        #         break

        final_el = elem[final_index]

        for i in range(len(final_el)):
            if final_el[i] == 'brBH':
                #print sites[i][0], sites[i][1]
                site_x = int(sites[i][0]*2.0)
                site_y = int(sites[i][1])
                #print site_x, site_y
                inserted_loc[dimer_window][site_x,site_y] += 1
                incorp_sites.append((sites[i][0]*0.7666,sites[i][1]*0.7666))

        for i in range(len(incorp_sites)):
            site_a = incorp_sites[i]
            for j in range(len(incorp_sites)):
                if i != j:
                    site_b = incorp_sites[j]
                    distance = math.sqrt((site_a[0]-site_b[0])**2 + (site_a[1]-site_b[1])**2)
                    distances.append(distance)

        incorp = ['u','BH','u']
        final_count = 0
        incorp_string = str(incorp).strip('[]')
        incorp_string = "'brBH'"
        second_incorp_string = "'BH'"
        #print incorp_string
        final_sites_string = str(elem[final_index]).strip('[]')
        #print final_sites_string
        inserted[dimer_window].append(final_sites_string.count(incorp_string)+final_sites_string.count(second_incorp_string))
        #inserted[dimer_window].append(final_count)
        #inserted[dimer_window].append(float(final_sites_string.count(incorp_string))/float(final_sites_string.count(total_p_string)))
    end_time = time.time()
    wall_time[dimer_window]=end_time-start_time

data = {}
data['map'] = inserted_loc[dimer_window].tolist()
data['nruns'] = nruns
data['distances'] = distances

json.dump(data,open('%sx%sincorp_map.json'%(dimer_window,dimer_window),'w'))
print 'inserted count', inserted
print 'wall times', wall_time

incorp_percent = {}

for dimer_window in patch_size:
    incorp_percent[dimer_window] = np.average(np.array(inserted[dimer_window])/num_avail_dimers(dimer_window,dimer_window))
    plt.imshow(np.matrix.transpose(inserted_loc[dimer_window]/nruns),interpolation='nearest',origin = 'lower',
                extent=(-2.*dimer_width,dimer_window+2.*dimer_width,-1.*dimer_height,dimer_window+1.*dimer_height))
    cbar = plt.colorbar()
    cbar.set_label('Proportion of B incorporation')
    plt.xlabel('(nm)')
    plt.ylabel('(nm)')
    plt.savefig('%s-%s_heatmap_anneal.pdf'%(dimer_window,dimer_window))
    #plt.show()
    plt.close()
    # planar_avg_incorp= np.sum(inserted_loc[dimer_window]/np.sum(inserted[dimer_window]),axis=1)
    # #print len(planar_avg_incorp)
    # #print len(np.arange(0.0,(len(planar_avg_incorp))*dimer_width,dimer_width))
    # plt.plot(np.arange(0.0,(len(planar_avg_incorp))*dimer_width,dimer_width),planar_avg_incorp)
    # plt.xlabel('X axis (nm)')
    # plt.ylabel('Proportion of P incorporation')
    # plt.savefig('%s-%s-planar_x_avg_anneal.jpg'%(dimer_window,dimer_window))
    # plt.close()
    #
    # planar_avg_incorp= np.sum(inserted_loc[dimer_window]/np.sum(inserted[dimer_window]),axis=0)
    # plt.plot(np.arange(0.0,len(planar_avg_incorp)*dimer_height,dimer_height),planar_avg_incorp)
    # plt.xlabel('Y axis (nm)')
    # plt.ylabel('Proportion of P incorporation')
    # plt.savefig('%s-%s-planar_y_avg_anneal.jpg'%(dimer_window,dimer_window))
    # plt.close()


print 'incorporated percent', incorp_percent
