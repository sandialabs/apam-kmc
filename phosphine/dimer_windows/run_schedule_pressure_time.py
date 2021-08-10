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

times = [6,60,600,6000,60000]
pressures = ['3e-12','3e-11','3e-10','3e-9','3e-8']

for time in times:
    for pressure in pressures:

        #for attempt_freq in [1e12,1e13,1e14,1e15,1e16]:#1e12,1e13,1e14,
        for attempt_freq in [1e12]:
            for dimer_window in [2,2.5,3,3.5,4]:#range(2,5):
            #for dimer_window in range(2,5):
            #for dimer_window in [3.5,4.5]:
            #for dimer_window in range(3,4):
                inserted[dimer_window] = []
                inserted_first[dimer_window] = []

                for i in range(nruns):
                    num_incorporated = 0
                    interactions  = KMCInteractionsFromScript('./rxns-pressure/schultz-reactions-A%s-P%s.py'%(attempt_freq,pressure))
                    control_parameters = KMCControlParameters(number_of_steps=int(3*10*dimer_window),
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


                    incorp = ['PH']
                    incorp_string = str(incorp).strip('[]')
                    final_sites_string = str(elem[final_index]).strip('[]')
                    inserted_first[dimer_window].append(final_sites_string.count(incorp_string))

                    print elem[final_index]
                    new_configuration = KMCConfigurationFromScript('./configs-time/dose%s/%sdimconfig_highT.py'%(time,dimer_window))
                    #new_configuration.__types = elem[final_index]

                    #print new_configuration.__sites

                    #print configuration
                    interactions  = KMCInteractionsFromScript('schultz-reactions-A%s-highT.py'%(attempt_freq))
                    control_parameters = KMCControlParameters(number_of_steps=int(80*dimer_window),
                                                                dump_interval=1,seed=random.randint(1,2e9))
                    new_model = KMCLatticeModel(new_configuration, interactions)
                    try:
                        new_model.run(control_parameters,trajectory_filename='%sdimwindow-T.py'%(dimer_window))
                    except:
                        print ( ' No more options available')

                    if os.path.exists('%sdimwindow-T.py'%(dimer_window)):
                        trajfile='%sdimwindow-T.py'%(dimer_window)
                        time_to_check_after = 15
                    else:
                        trajfile="%sdimwindow.py"%(dimer_window)
                        time_to_check_after = time
                    global_dict = {}
                    local_dict  = {}
                    execfile(trajfile, global_dict, local_dict)
                    times = local_dict['times']
                    sites = local_dict["sites"]
                    elem = local_dict["types"]


                    #if os.path.exists('%sdimwindow-T.py'%(dimer_window)):
                    #    os.remove('%sdimwindow-T.py'%(dimer_window))


                    final_index = -1
                    for j in range(len(times)):
                        if times[j] > time_to_check_after:
                            final_index = j
                            break


                    incorp = ['PH']
                    incorp_string = str(incorp).strip('[]')
                    final_sites_string = str(elem[final_index]).strip('[]')
                    inserted[dimer_window].append(final_sites_string.count(incorp_string))

                    if inserted_first[dimer_window][-1] == 1 and inserted[dimer_window][-1] == 0:
                        break

                    # #print (inserted[dimer_window])
                    # if inserted[dimer_window][-1] == 0:
                    #     break


                print 'initial incorporated rate ', inserted_first
                print 'incorporated rate', inserted


            incorporated = {}
            g = open('incorporation-A%s-%s-runs-pressure%s-time%s.csv'%(attempt_freq,nruns,pressure,time),'w')
            g.write("Dimer Window,     0 inserted %%,    1 inserted %%,  2 inserted %%, 3 inserted %% \n")


            #for i in range(2,5):
            for i in [2,2.5,3,3.5,4]:
                incorporated[i] = [inserted[i].count(0)/float(nruns),inserted[i].count(1)/float(nruns),
                                    inserted[i].count(2)/float(nruns),inserted[i].count(3)/float(nruns)]
                g.write('%s,  %s,  %s,  %s,  %s \n'%(i,incorporated[i][0],incorporated[i][1],incorporated[i][2],incorporated[i][3]))

            g.close()
            # ninsert = np.arange(4)
            # width = 0.2
            # ymax = 0.8
            # title = "KMClib Insertion Statistics"
            #
            # fig,ax = plt.subplots(figsize=(12,10))
            # ax.bar(ninsert-width,incorporated[2],width,color='b',label='2 dimers')
            # ax.bar(ninsert,incorporated[3],width,color='r',label='3 dimers')
            # ax.bar(ninsert+width,incorporated[4],width,color='lightblue',label='4 dimers')
            # #ax.bar(ninsert+2*width,incorporated[5],width,color='b',label='5 dimers')
            # #ax.bar(ninsert+3*width,incorporated[6],width,color='g',label='6 dimers')
            # ax.set_xticks(ninsert+2*width)
            # ax.set_xticklabels(ninsert)
            # ax.legend(bbox_to_anchor=(1, 1), borderaxespad=0.)
            # ax.set_ylabel('fractional occurance')
            # ax.set_xlabel('number of features per site')
            # ax.set_title(title)
            # ax.axis(ymax=ymax,xmin=-0.2)
            # plt.savefig('kmclib_insertion_histogram_schultz-A%s-%s-h-desorption.png'%(attempt_freq,nruns))
            # plt.show()
