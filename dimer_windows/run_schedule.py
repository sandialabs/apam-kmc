#Quinn Campbell qcampbe@sandia.gov

from KMCLib import *
import matplotlib.pyplot as plt
import random
import numpy as np
import matplotlib as mpl

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Avenir'

nruns = 1

inserted= {}
for dimer_window in range(3,7):
#for dimer_window in range(3,4):
    inserted[dimer_window] = []
    for i in range(nruns):
        num_incorporated = 0
        interactions  = KMCInteractionsFromScript('diborane-rxn.py')
        control_parameters = KMCControlParameters(number_of_steps=30,
                                                    dump_interval=1,seed=random.randint(1,2e9))
        configuration = KMCConfigurationFromScript('%sdimconfig.py'%dimer_window)
        model = KMCLatticeModel(configuration,interactions)
        try:
           model.run(control_parameters,trajectory_filename='%sdimwindow.py'%(dimer_window))
        except:
            print ( ' No more options available')


        trajfile="%sdimwindow.py"%(dimer_window)
        global_dict = {}
        local_dict  = {}
        execfile(trajfile, global_dict, local_dict)
        times = local_dict['times']
        sites = local_dict["sites"]
        elem = local_dict["types"]

        final_index = -1
        time_to_check_after = 600.0  #In units of s
        for j in range(len(times)):
            if times[j] > time_to_check_after:
                final_index = j
                break

        incorp = ['u','BH','u']
        incorp_string = str(incorp).strip('[]')
        incorp_string = "'brBH'"
        #print incorp_string
        final_sites_string = str(elem[final_index]).strip('[]')
        #print final_sites_string
        #inserted[dimer_window].append(final_sites_string.count(incorp_string))
        # if inserted[dimer_window][-1] == 2:
        #     break


        new_configuration = KMCConfigurationFromScript('%sdimconfig_highT.py'%dimer_window)
        #new_configuration.__types = elem[final_index]

        #print new_configuration.__sites

        #print configuration
        interactions  = KMCInteractionsFromScript('diborane-rxn-highT.py')
        control_parameters = KMCControlParameters(number_of_steps=100,#*dimer_window,
                                                    dump_interval=1,seed=random.randint(1,2e9))
        new_model = KMCLatticeModel(new_configuration, interactions)
        try:
            new_model.run(control_parameters,trajectory_filename='%sdimwindow-T.py'%(dimer_window))
        except:
            print ( ' No more options available')

        trajfile="%sdimwindow-T.py"%(dimer_window)
        global_dict = {}
        local_dict  = {}
        execfile(trajfile, global_dict, local_dict)
        times = local_dict['times']
        sites = local_dict["sites"]
        elem = local_dict["types"]

        final_index = -1
        time_to_check_after = 60.0  #In units of s
        for j in range(len(times)):
            if times[j] > time_to_check_after:
                final_index = j
                break

        incorp = ['u','BH','u']
        incorp_string = str(incorp).strip('[]')
        incorp_string = "'brBH'"
        #print incorp_string
        final_sites_string = str(elem[final_index]).strip('[]')
        #print final_sites_string
        inserted[dimer_window].append(final_sites_string.count(incorp_string))

    print 'incorporated rate', inserted


incorporated = {}
g = open('incorporation-%s-runs-rxn-down-test.csv'%nruns,'w')
g.write("Dimer Window,     0 inserted %%,    1 inserted %%,  2 inserted %%, 3 inserted %% \n")


for i in range(3,7):
    incorporated[i] = [inserted[i].count(0)/float(nruns),inserted[i].count(1)/float(nruns),
                        inserted[i].count(2)/float(nruns),inserted[i].count(3)/float(nruns),inserted[i].count(4)/float(nruns),
                        inserted[i].count(5)/float(nruns),inserted[i].count(6)/float(nruns)]
    g.write('%s,  %s,  %s,  %s,  %s,  %s,  %s,  %s \n'%(i,incorporated[i][0],incorporated[i][1],incorporated[i][2],incorporated[i][3],incorporated[i][4],incorporated[i][5],incorporated[i][6]))

ninsert = np.arange(7)
width = 0.2
ymax = 0.8
title = "KMClib Insertion Statistics"

fig,ax = plt.subplots(figsize=(12,10))
ax.bar(ninsert,incorporated[3],width,color='r',label='3 dimers')
ax.bar(ninsert+width,incorporated[4],width,color='lightblue',label='4 dimers')
ax.bar(ninsert+2*width,incorporated[5],width,color='b',label='5 dimers')
ax.bar(ninsert+3*width,incorporated[6],width,color='g',label='6 dimers')
#ax.bar(ninsert+4*width,incorporated[7],width,color='m',label='7 dimers')
ax.set_xticks(ninsert+2*width)
ax.set_xticklabels(ninsert)
ax.legend(bbox_to_anchor=(1, 1), borderaxespad=0.)
ax.set_ylabel('fractional occurance')
ax.set_xlabel('number of features per site')
ax.set_title(title)
ax.axis(ymax=ymax,xmin=-0.2)
plt.savefig('diborane-%s-sched-desorption.png'%nruns)
plt.show()
