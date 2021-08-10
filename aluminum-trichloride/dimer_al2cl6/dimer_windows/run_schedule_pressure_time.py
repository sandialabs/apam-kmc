from KMCLib import *
import matplotlib.pyplot as plt
import random
import numpy as np
import matplotlib as mpl
import os

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Avenir'

nruns =200

inserted= {}
inserted_first = {}

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
                time_to_check_after = 900  #In units of s
                for j in range(len(times)):
                    if times[j] > time_to_check_after:
                        final_index = j-1
                        break


                incorp = ['brAlCl']
                incorp_string = str(incorp).strip('[]')
                final_sites_string = str(elem[final_index]).strip('[]')
                inserted_first[dimer_window].append(final_sites_string.count(incorp_string)/2)

                print elem[final_index]
                new_configuration = KMCConfigurationFromScript('./configs-time/%sdimconfig-highT-t%s.py'%(dimer_window,time))
                #new_configuration.__types = elem[final_index]

                #print new_configuration.__sites

                #print configuration
                interactions  = KMCInteractionsFromScript('alcl3-rxn-br-modify-high-T.py')
                control_parameters = KMCControlParameters(number_of_steps=10*dimer_window,
                                                            dump_interval=1,seed=random.randint(1,2e9))
                new_model = KMCLatticeModel(new_configuration, interactions)
                try:
                    new_model.run(control_parameters,trajectory_filename='%sdimwindow-T.py'%(dimer_window))
                except:
                    print ( ' No more options available')

                if os.path.exists('%sdimwindow-T.py'%(dimer_window)):
                    trajfile='%sdimwindow-T.py'%(dimer_window)
                    time_to_check_after = 5
                else:
                    trajfile="%sdimwindow.py"%(dimer_window)
                    time_to_check_after = time
                global_dict = {}
                local_dict  = {}
                execfile(trajfile, global_dict, local_dict)
                times = local_dict['times']
                sites = local_dict["sites"]
                elem = local_dict["types"]


                if os.path.exists('%sdimwindow-T.py'%(dimer_window)):
                   os.remove('%sdimwindow-T.py'%(dimer_window))


                final_index = -1
                for j in range(len(times)):
                    if times[j] > time_to_check_after:
                        final_index = j
                        break


                incorp = ['brAlCl']
                incorp_string = str(incorp).strip('[]')
                final_sites_string = str(elem[final_index]).strip('[]')
                #print (final_sites_string)
                dimer_incorp = ['brAlCl','u','brAlCl','brAlCl','u','brAlCl']
                dimer_incorp_string = str(dimer_incorp).strip('[]')
                strict_incorp1 = ['brAlCl','u','Cl','brAlCl','u','AlCl2']
                strict_incorp1_string = str(strict_incorp1).strip('[]')
                #print (strict_incorp1_string)
                strict_incorp2 = ['brAlCl','u','AlCl2','brAlCl','u','Cl']
                strict_incorp2_string = str(strict_incorp2).strip('[]')
                strict_incorp3 = ['Cl','u','brAlCl','AlCl2','u','brAlCl']
                strict_incorp3_string = str(strict_incorp3).strip('[]')
                strict_incorp4 = ['AlCl2','u','brAlCl','Cl','u','brAlCl']
                strict_incorp4_string = str(strict_incorp4).strip('[]')
                total_subtract = 2*final_sites_string.count(dimer_incorp_string) + final_sites_string.count(strict_incorp1_string)+ final_sites_string.count(strict_incorp2_string)+ final_sites_string.count(strict_incorp3_string)+ final_sites_string.count(strict_incorp4_string)
                inserted[dimer_window].append(final_sites_string.count(incorp_string)/2 - total_subtract)

                # if inserted_first[dimer_window][-1] == 1 and inserted[dimer_window][-1] == 0:
                #     break

                # #print (inserted[dimer_window])
                # if inserted[dimer_window][-1] == 0:
                #     break


            print 'initial incorporated rate ', inserted_first
            print 'incorporated rate', inserted


        incorporated = {}
        g = open('incorporation-%s-runs-single-acceptor-P-%s-t-%s-strict-count.csv'%(nruns,pressure,time),'w')
        g.write("Dimer Window,     0 inserted %%,    1 inserted %%,  2 inserted %% \n")


        for i in range(3,4):
            incorporated[i] = [inserted[i].count(0)/float(nruns),inserted[i].count(1)/float(nruns),
                                inserted[i].count(2)/float(nruns)]
            g.write('%s,  %s,  %s,  %s \n'%(i,incorporated[i][0],incorporated[i][1],incorporated[i][2]))
#
# ninsert = np.arange(4)
# width = 0.2
# ymax = 0.8
# title = "KMClib Insertion Statistics"
#
# fig,ax = plt.subplots(figsize=(12,10))
# ax.bar(ninsert,incorporated[3],width,color='r',label='3 dimers')
# ax.bar(ninsert+width,incorporated[4],width,color='lightblue',label='4 dimers')
# ax.bar(ninsert+2*width,incorporated[5],width,color='b',label='5 dimers')
# #ax.bar(ninsert+3*width,incorporated[6],width,color='g',label='6 dimers')
# ax.set_xticks(ninsert+2*width)
# ax.set_xticklabels(ninsert)
# ax.legend(bbox_to_anchor=(1, 1), borderaxespad=0.)
# ax.set_ylabel('fractional occurance')
# ax.set_xlabel('number of features per site')
# ax.set_title(title)
# ax.axis(ymax=ymax,xmin=-0.2)
# plt.savefig('kmclib_insertion_histogram_schultz-%s-h-desorption.png'%(nruns))
# plt.show()
