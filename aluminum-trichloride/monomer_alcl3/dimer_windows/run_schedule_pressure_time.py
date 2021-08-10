from KMCLib import *
import matplotlib.pyplot as plt
import random
import numpy as np
import matplotlib as mpl
import os

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Avenir'

nruns = 200

inserted= {}
inserted_first = {}

#for dimer_window in range(2,4):
times = [9,90,900,9000,90000]
pressures = ['4e-11','4e-10','4e-9','4e-8']

for time in times:
    for pressure in pressures:
        for dimer_window in range(3,4):
            inserted[dimer_window] = []
            inserted_first[dimer_window] = []

            for i in range(nruns):
                num_incorporated = 0
                interactions  = KMCInteractionsFromScript('./rxns-pressure/alcl3-rxn-br-modify-P%s.py'%pressure)
                control_parameters = KMCControlParameters(number_of_steps=10*dimer_window,
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
                time_to_check_after = time  #In units of s
                for j in range(len(times)):
                    if times[j] > time_to_check_after:
                        final_index = j
                        break


                incorp = ['brAlCl']
                incorp_string = str(incorp).strip('[]')
                final_sites_string = str(elem[final_index]).strip('[]')
                inserted_first[dimer_window].append(final_sites_string.count(incorp_string)/2)

                print elem[final_index]
                # new_configuration = KMCConfigurationFromScript('%sdimconfig_highT.py'%dimer_window)
                # #new_configuration.__types = elem[final_index]
                #
                # #print new_configuration.__sites
                #
                # #print configuration
                # interactions  = KMCInteractionsFromScript('alcl3-rxn-br-modify-high-T.py')
                # control_parameters = KMCControlParameters(number_of_steps=20*dimer_window,
                #                                             dump_interval=1,seed=random.randint(1,2e9))
                # new_model = KMCLatticeModel(new_configuration, interactions)
                # try:
                #     new_model.run(control_parameters,trajectory_filename='%sdimwindow-T.py'%(dimer_window))
                # except:
                #     print ( ' No more options available')
                #
                # if os.path.exists('%sdimwindow-T.py'%(dimer_window)):
                #     trajfile='%sdimwindow-T.py'%(dimer_window)
                #     time_to_check_after = 5
                # else:
                #     trajfile="%sdimwindow.py"%(dimer_window)
                #     time_to_check_after = 900
                # global_dict = {}
                # local_dict  = {}
                # execfile(trajfile, global_dict, local_dict)
                # times = local_dict['times']
                # sites = local_dict["sites"]
                # elem = local_dict["types"]
                #
                #
                # if os.path.exists('%sdimwindow-T.py'%(dimer_window)):
                #    os.remove('%sdimwindow-T.py'%(dimer_window))


                # final_index = -1
                # for j in range(len(times)):
                #     if times[j] > time_to_check_after:
                #         final_index = j-1
                #         break


                incorp = ['brAlCl']
                incorp_string = str(incorp).strip('[]')
                final_sites_string = str(elem[final_index]).strip('[]')
                inserted[dimer_window].append(final_sites_string.count(incorp_string)/2)

                if inserted_first[dimer_window][-1] == 1 and inserted[dimer_window][-1] == 0:
                    break

                # #print (inserted[dimer_window])
                # if inserted[dimer_window][-1] == 0:
                #     break


            print 'initial incorporated rate ', inserted_first
            print 'incorporated rate', inserted


        incorporated = {}
        g = open('incorporation-%s-runs-single-acceptor-P-%s-t-%s.csv'%(nruns,pressure,time),'w')
        g.write("Dimer Window,     0 inserted %%,    1 inserted %%,  2 inserted %%\n")


        for i in range(3,4):
            incorporated[i] = [inserted[i].count(0)/float(nruns),inserted[i].count(1)/float(nruns),
                                inserted[i].count(2)/float(nruns)]
            g.write('%s,  %s,  %s,  %s \n'%(i,incorporated[i][0],incorporated[i][1],incorporated[i][2]))
