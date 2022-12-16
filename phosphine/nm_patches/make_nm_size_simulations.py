import os
import subprocess
import shutil

def main():

    sizes = [0.5,1.0,2.0,3.0,5.0,10.0,15.0,20.0]


    for endpt in [,'brPH','allP']:
        os.mkdir(endpt)
        os.chdir(endpt)
        for size in sizes:
            os.mkdir('{}x{}nm'.format(size,size))
            os.chdir('{}x{}nm'.format(size,size))

            make_config_files(size)
            make_schedule_file(size,endpt)
            make_submit_file(size,endpt)

            if endpt == 'allP':
                shutil.copy('../../cross-reactions-A1e+12-PH.py','./cross-reactions-A1e+12.py')
                shutil.copy('../../cross-reactions-A1e+12-highT-PH.py','./cross-reactions-A1e+12-highT.py')
            else:
                shutil.copy('../../cross-reactions-A1e+12-P.py','./cross-reactions-A1e+12.py')
                shutil.copy('../../cross-reactions-A1e+12-highT-P.py','./cross-reactions-A1e+12-highT.py')

            os.mkdir('incorp_kmc_data')
            os.mkdir('heatmaps')

            subprocess.call(['sbatch','submit.bash'])

            os.chdir('..')

        os.chdir('..')

    return

def make_config_files(size):

    f = open('3-3config.py','w')

    f.write("""#
# This is the lattice setup for the Schultz mapping.
# Key finding: PH3 only likely to adsorb in up on one dimer, down on
# the next dimer fashion. This geometry seems to be the main characteristic
# change. Requires a lattice that takes into account two dimer sites.
#

from KMCLib import *
import math

# ----------------------------------------------------------------------------
# Unit cell

cell_vectors = [[   2.000000e+00,   0.000000e+00,   0.000000e+00],
                [   0.000000e+00,   6.000000e+00,   0.000000e+00],
                [   0.000000e+00,   0.000000e+00,   1.000000e+00]]

#Here 0.00 refers to the first dimer site
basis_points = [[   0.000000e+00,   0.000000e+00/2.0,   0.000000e+00],
               [   0.000000e+00,   0.333333e+00/2.0,   0.000000e+00],
               [   0.000000e+00,   0.666666e+00/2.0,   0.000000e+00],
               [   0.500000e+00,   0.000000e+00/2.0,   0.000000e+00],
               [   0.500000e+00,   0.333333e+00/2.0,   0.000000e+00],
               [   0.500000e+00,   0.666666e+00/2.0,   0.000000e+00], #all the sites after this will be a row up
               [   0.000000e+00,   0.000000e+00/2.0+0.5,   0.000000e+00],
               [   0.000000e+00,   0.333333e+00/2.0+0.5,   0.000000e+00],
               [   0.000000e+00,   0.666666e+00/2.0+0.5,   0.000000e+00],
               [   0.500000e+00,   0.000000e+00/2.0+0.5,   0.000000e+00],
               [   0.500000e+00,   0.333333e+00/2.0+0.5,   0.000000e+00],
               [   0.500000e+00,   0.666666e+00/2.0+0.5,   0.000000e+00]]



unit_cell = KMCUnitCell(
    cell_vectors=cell_vectors,
    basis_points=basis_points)

# -----------------------------------------------------------------------------
# Lattice

nm_width = {}
nm_height = {} 

dimer_height = 0.766
dimer_width = 0.3874

num_dimers_wide = int(math.ceil(nm_width/dimer_width))
num_dimers_tall = int(math.ceil(nm_height/dimer_height))


reps_x = int(math.ceil((num_dimers_wide+3.)/2.))
reps_y = int(math.ceil((num_dimers_tall)/2.) + 1)
print reps_y

lattice = KMCLattice(
    unit_cell=unit_cell,
    repetitions=(reps_x,reps_y,1),
    periodic=(False, False, False))


#print lattice.sites()
types = []
#Having a for loop seems like it will be simpler  than trying to figure out how sites
# reorder on my own every single time (presuming you don't mind long or staments)
for site in lattice.sites():
    #print site
    if (site[0] < 1.0) or (site[0] >= (num_dimers_wide +2)*0.5) or (site[1]<1.0) or (site[1] >= (num_dimers_tall+2.)*0.5):
        ending_decimal = site[1]-math.floor(site[1])
        #print ending_decimal
        if abs(ending_decimal - 0.333333/2.0) < 0.05 or abs(ending_decimal - 0.333333e+00/2.0-0.5) < 0.05:
            #print 'u'
            types.append('u')
        else:
            #print 'H'
            types.append('H')
    else:
        types.append('u')


#print types

#I will use 'u' for unoccupied
# and 'b' for involved in a bridging reaction
# otherwise types will match their atomic label (i.e PH2)
#I'm including a few types such as H2 and P that we essentially won't run into
#Shouldn't be an issue, but if it's causing issues later, here's your problem
possible_types = ['u','b','PH3','PH2','PH','P','H','H2']

configuration = KMCConfiguration(
    lattice=lattice,
    types=types,
    possible_types=possible_types)
""".format(size,size))

    f.close()

    g = open('3-3config-highT.py','w')

    g.write("""#
# This is the lattice setup for the Schultz mapping.
# Key finding: PH3 only likely to adsorb in up on one dimer, down on
# the next dimer fashion. This geometry seems to be the main characteristic
# change. Requires a lattice that takes into account two dimer sites.
#

from KMCLib import *
import pickle
import math
import json
import numpy as np

# ----------------------------------------------------------------------------
# Unit cell

cell_vectors = [[   2.00e+00,   0.000000e+00,   0.000000e+00],
                [   0.000000e+00,   3.0000e+00,   0.000000e+00],
                [   0.000000e+00,   0.000000e+00,   1.000000e+00]]

#Here 0.00 refers to the first dimer site
cell_vectors = [[   2.000000e+00,   0.000000e+00,   0.000000e+00],
                [   0.000000e+00,   6.000000e+00,   0.000000e+00],
                [   0.000000e+00,   0.000000e+00,   1.000000e+00]]

#Here 0.00 refers to the first dimer site
basis_points = [[   0.000000e+00,   0.000000e+00/2.0,   0.000000e+00],
               [   0.000000e+00,   0.333333e+00/2.0,   0.000000e+00],
               [   0.000000e+00,   0.666666e+00/2.0,   0.000000e+00],
               [   0.500000e+00,   0.000000e+00/2.0,   0.000000e+00],
               [   0.500000e+00,   0.333333e+00/2.0,   0.000000e+00],
               [   0.500000e+00,   0.666666e+00/2.0,   0.000000e+00], #all the sites after this will be a row up
               [   0.000000e+00,   0.000000e+00/2.0+0.5,   0.000000e+00],
               [   0.000000e+00,   0.333333e+00/2.0+0.5,   0.000000e+00],
               [   0.000000e+00,   0.666666e+00/2.0+0.5,   0.000000e+00],
               [   0.500000e+00,   0.000000e+00/2.0+0.5,   0.000000e+00],
               [   0.500000e+00,   0.333333e+00/2.0+0.5,   0.000000e+00],
               [   0.500000e+00,   0.666666e+00/2.0+0.5,   0.000000e+00]]



unit_cell = KMCUnitCell(
    cell_vectors=cell_vectors,
    basis_points=basis_points)

# -----------------------------------------------------------------------------
# Lattice
nm_width = {}
nm_height = {}

dimer_height = 0.766
dimer_width = 0.3874

num_dimers_wide = int(math.ceil(nm_width/dimer_width))
num_dimers_tall = int(math.ceil(nm_height/dimer_height))


reps_x = int(math.ceil((num_dimers_wide+3.)/2.))
reps_y = int(math.ceil((num_dimers_tall)/2.) + 1)
#print reps

lattice = KMCLattice(
    unit_cell=unit_cell,
    repetitions=(reps_x,reps_y,1),
    periodic=(False, False, False))


#reading in from the earlier incorp
trajfile="./dimwindow.py"
global_dict = {{}}
local_dict  = {{}}
execfile(trajfile, global_dict, local_dict)
times = local_dict['times']
sites = local_dict["sites"]
elem = local_dict["types"]

final_index = -1
time_to_check_after = 480.0  #In units of s
for j in range(len(times)):
    if times[j] > time_to_check_after:
        final_index = j-1
        print final_index
        break


types = elem[final_index]


##types = ['u']=
#print types

#I will use 'u' for unoccupied
# and 'b' for involved in a bridging reaction
# otherwise types will match their atomic label (i.e PH2)
#I'm including a few types such as H2 and P that we essentially won't run into
#Shouldn't be an issue, but if it's causing issues later, here's your problem
possible_types = ['u','b','PH3','PH2','PH','P','H','H2']

configuration = KMCConfiguration(
    lattice=lattice,
    types=types,
    possible_types=possible_types)""".format(size,size))

    g.close()

    return

def make_schedule_file(size,endpt):

    if endpt == 'allP':
        f = open('run_schedule.py','w')

        f.write("""from KMCLib import *
import matplotlib.pyplot as plt
import random
import numpy as np
import matplotlib as mpl
import os
import json
import math

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Arial'

nruns = 200
total_nruns=200

inserted = []
inserted_first = []
nm_width = {}
nm_height = {}

dimer_height = 0.766
dimer_width = 0.3874

num_dimers_wide = int(math.ceil(nm_width/dimer_width))
num_dimers_tall = int(math.ceil(nm_height/dimer_height))

inserted_loc = np.zeros((num_dimers_wide+3,num_dimers_tall+3))

data = {{}}
data['runs'] =  0

for i in range(nruns):
    if os.path.isfile('./incorp_kmc_data/incorp.json'):
        data = json.load(open('./incorp_kmc_data/incorp.json','r'))

        inserted_loc = np.array(data['map'])
        inserted = data['inserted']

    if data['runs'] > total_nruns:
        break

    num_incorporated = 0
    data['runs'] +=1
    print data['runs'], " runs so far"
    interactions  = KMCInteractionsFromScript('cross-reactions-A1e+12.py')
    control_parameters = KMCControlParameters(number_of_steps=int(100*nm_width*nm_height),
                                                dump_interval=10,seed=random.randint(1,2e9))
    configuration = KMCConfigurationFromScript('./3-3config.py')
    model = KMCLatticeModel(configuration,interactions)
    try:
        model.run(control_parameters,trajectory_filename='dimwindow.py')
    except:
        print ( ' No more options available')


    trajfile="dimwindow.py"
    global_dict = {{}}
    local_dict  = {{}}
    execfile(trajfile, global_dict, local_dict)
    times = local_dict['times']
    sites = local_dict["sites"]
    elem = local_dict["types"]

    final_index = -1
    time_to_check_after = 480.0  #In units of s
    for j in range(len(times)):
        if times[j] > time_to_check_after:
            final_index = j-1
            break


    incorp = 'P'
    incorp_string = str(incorp)
    final_sites_string = str(elem[final_index]).strip('[]')
    inserted_first.append(final_sites_string.count(incorp_string))

    print elem[final_index]
    new_configuration = KMCConfigurationFromScript('3-3config-highT.py')
    #new_configuration.__types = elem[final_index]

    #print new_configuration.__sites

    #print configuration
    interactions  = KMCInteractionsFromScript('cross-reactions-A1e+12-highT.py')
    control_parameters = KMCControlParameters(number_of_steps=int(2000*nm_width*nm_height),
                                                dump_interval=100,seed=random.randint(1,2e9))
    new_model = KMCLatticeModel(new_configuration, interactions)
    try:
        new_model.run(control_parameters,trajectory_filename='dimwindow-T.py')
    except:
        print ( ' No more options available')

    if os.path.exists('dimwindow-T.py'):
        trajfile='dimwindow-T.py'
        time_to_check_after = 10
    else:
        trajfile="dimwindow.py"
        time_to_check_after = 480
    global_dict = {{}}
    local_dict  = {{}}
    execfile(trajfile, global_dict, local_dict)
    times = local_dict['times']
    sites = local_dict["sites"]
    elem = local_dict["types"]


    #if os.path.exists('%sdimwindow-T.py'%(dimer_window)):
    #    os.remove('%sdimwindow-T.py'%(dimer_window))


    final_index = -1
    for j in range(len(times)):
        if times[j] > time_to_check_after:
            final_index = j-1
            break
    print final_index


    final_el = elem[final_index]
    incorp = 'P'
    incorp_string = str(incorp)
    final_sites_string = str(final_el).strip('[]')
    inserted.append(final_sites_string.count(incorp_string))


    print len(sites)
    for i in range(len(final_el)):
        if final_el[i] == 'PH':
            print sites[i][0], sites[i][1]
            site_x = int(sites[i][0]*2.0)
            site_y = (int(sites[i][1]*2.0))#+2
            #print site_x, site_y
            inserted_loc[site_x,site_y] += 1
            #inserted_loc[site_x,site_y+1] += 1

    data['map'] = inserted_loc.tolist()
    data['inserted'] = inserted
    # if inserted[-1] < 3:
    #     break

    json.dump(data,open('./incorp_kmc_data/incorp.json','w'))

    # #print (inserted[dimer_window])
    # if inserted[dimer_window][-1] == 0:
    #     break


print 'initial incorporated rate ', inserted_first
print 'incorporated rate', inserted
print data['runs']

plt.imshow(np.transpose(inserted_loc/float(data['runs'])),interpolation='nearest',origin = 'lower',aspect=dimer_height/dimer_width)#,
            #extent=(0,dimer_window+3.*dimer_width,0,dimer_window+3.*dimer_height))
cbar = plt.colorbar()
cbar.set_label('Proportion of P incorporation')
plt.xlabel('(position)')
plt.ylabel('(position)')
plt.savefig('./heatmaps/heatmap_anneal.jpg')
#plt.show()
plt.close()""".format(size,size))

        f.close()
    else:
        f = open('run_schedule.py','w')

        f.write("""from KMCLib import *
import matplotlib.pyplot as plt
import random
import numpy as np
import matplotlib as mpl
import os
import json
import math

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Arial'

nruns = 500
total_nruns=500

inserted = []
inserted_first = []
nm_width = {}
nm_height = {}

dimer_height = 0.766
dimer_width = 0.3874

num_dimers_wide = int(math.ceil(nm_width/dimer_width))
num_dimers_tall = int(math.ceil(nm_height/dimer_height))

inserted_loc = np.zeros((num_dimers_wide+3,num_dimers_tall+3))

data = {{}}
data['runs'] =  0

num_not_finished = 0

for i in range(nruns):
    if os.path.isfile('./incorp_kmc_data/incorp.json'):
        data = json.load(open('./incorp_kmc_data/incorp.json','r'))

        inserted_loc = np.array(data['map'])
        inserted = data['inserted']

    if data['runs'] > total_nruns:
        break

    num_incorporated = 0
    data['runs'] +=1
    print data['runs'], " runs so far"
    interactions  = KMCInteractionsFromScript('cross-reactions-A1e+12.py')
    control_parameters = KMCControlParameters(number_of_steps=int(20*nm_height*nm_width),
                                                dump_interval=1,seed=random.randint(1,2e9))
    configuration = KMCConfigurationFromScript('./3-3config.py')
    model = KMCLatticeModel(configuration,interactions)
    try:
        model.run(control_parameters,trajectory_filename='dimwindow.py')
    except:
        print ( ' No more options available')


    trajfile="dimwindow.py"
    global_dict = {{}}
    local_dict  = {{}}
    execfile(trajfile, global_dict, local_dict)
    times = local_dict['times']
    sites = local_dict["sites"]
    elem = local_dict["types"]

    final_index = -1
    time_to_check_after = 480.0  #In units of s
    for j in range(len(times)):
        if times[j] > time_to_check_after:
            final_index = j-1
            break


    incorp = ['P']
    incorp_string = str(incorp).strip('[]')
    final_sites_string = str(elem[final_index]).strip('[]')
    inserted_first.append(final_sites_string.count(incorp_string))

    print elem[final_index]
    new_configuration = KMCConfigurationFromScript('3-3config-highT.py')
    #new_configuration.__types = elem[final_index]

    #print new_configuration.__sites

    #print configuration
    interactions  = KMCInteractionsFromScript('cross-reactions-A1e+12-highT.py')
    control_parameters = KMCControlParameters(number_of_steps=int(200*nm_height*nm_width),
                                                dump_interval=1,seed=random.randint(1,2e9))
    new_model = KMCLatticeModel(new_configuration, interactions)
    try:
        new_model.run(control_parameters,trajectory_filename='dimwindow-T.py')
    except:
        print ( ' No more options available')

    if os.path.exists('dimwindow-T.py'):
        trajfile='dimwindow-T.py'
        time_to_check_after = 10
    else:
        trajfile="dimwindow.py"
        time_to_check_after = 480
    global_dict = {{}}
    local_dict  = {{}}
    execfile(trajfile, global_dict, local_dict)
    times = local_dict['times']
    sites = local_dict["sites"]
    elem = local_dict["types"]


    #if os.path.exists('%sdimwindow-T.py'%(dimer_window)):
    #    os.remove('%sdimwindow-T.py'%(dimer_window))


    final_index = -1
    for j in range(len(times)):
        if times[j] > time_to_check_after:
            final_index = j-1
            break
    print final_index

    if final_index == -1:
        num_not_finished += 1


    final_el = elem[final_index]
    incorp = ['P']
    incorp_string = str(incorp).strip('[]')
    final_sites_string = str(final_el).strip('[]')
    inserted.append(final_sites_string.count(incorp_string))


    print len(sites)
    for i in range(len(final_el)):
        if final_el[i] == 'P':
            print sites[i][0], sites[i][1]
            site_x = int(sites[i][0]*2.0)
            site_y = (int(sites[i][1]*2.0))#+2
            #print site_x, site_y
            inserted_loc[site_x,site_y] += 1
            #inserted_loc[site_x,site_y+1] += 1

    data['map'] = inserted_loc.tolist()
    data['inserted'] = inserted
    # if inserted[-1] < 3:
    #     break

    json.dump(data,open('./incorp_kmc_data/incorp.json','w'))

    # #print (inserted[dimer_window])
    # if inserted[dimer_window][-1] == 0:
    #     break


print 'initial incorporated rate ', inserted_first
print 'incorporated rate', inserted
print data['runs']

print 'Number anneals not finished: ',num_not_finished

plt.imshow(np.transpose(inserted_loc/float(data['runs'])),interpolation='nearest',origin = 'lower',aspect=dimer_height/dimer_width)#,
            #extent=(0,dimer_window+3.*dimer_width,0,dimer_window+3.*dimer_height))
cbar = plt.colorbar()
cbar.set_label('Proportion of P incorporation')
plt.xlabel('(position)')
plt.ylabel('(position)')
plt.savefig('./heatmaps/heatmap_anneal.jpg')
#plt.show()
plt.close()
""".format(size,size))
        f.close()
    return

def make_submit_file(size,endpt):

    f = open('submit.bash','w')

    f.write("""#!/bin/bash
## Do not put any commands or blank lines before the #SBATCH lines
#SBATCH --nodes=1                     # Number of nodes - all cores per node are allocated to the job
#SBATCH --time=4:00:00                # Wall clock time (HH:MM:SS) - once the job exceeds this time, the job will be terminated (default is 5 minutes)
#SBATCH --account=FY190048            # WC ID
#SBATCH --job-name=P9{}{}nm               # Name of job
#SBATCH --partition=short,batch             # partition/queue name: short or batch
                                      #            short: 4hrs wallclock limit
                                      #            batch: nodes reserved for > 4hrs (default)
                                      #            short,batch: job will move to batch if short is full

#SBATCH --qos=normal                  # Quality of Service: long, large, priority or normal
                                      #           normal: request up to 48hrs wallclock (default)
                                      #           long:   request up to 96hrs wallclock and no larger than 64nodes
                                      #           large:  greater than 50% of cluster (special request)
                                      #           priority: High priority jobs (special request)

nodes=$SLURM_JOB_NUM_NODES           # Number of nodes - the number of nodes you have requested (for a list of SLURM environment variables see "man sbatch")
cores=1                             # Number MPI processes to run on each node (a.k.a. PPN)
                                     # tlcc2 has 16 cores per node
# using openmpi-intel/1.10
#module load mkl/16.0
#module load intel
#module load mkl
#module load intel-mpi
#module load fftw

/ascldap/users/qcampbe/anaconda2/bin/python run_schedule.py  | tee kmc.out 

# Note 1: This will start ($nodes * $cores) total MPI processes using $cores per node.
#           If you want some other number of processes, add "-np N" after the mpiexec, where N is the total you want.
#           Example:  mpiexec -np 24  ......(for a 2 node job, this will load 16 processes on the first node and 8 processes on the second node)
#           If you want a specific number of process to run on each node, (thus increasing the effective memory per core), use the --npernode option.
#           Example: mpiexec -np 24 --npernode 12  ......(for a 2 node job, this will load 12 processes on each node)

# The default version of Open MPI is version 1.10.

# For openmpi 1.10: mpiexec --bind-to core --npernode 8 --n PUT_THE_TOTAL_NUMBER_OF_MPI_PROCESSES_HERE /path/to/executable [--args...]

# To submit your job, do:
# sbatch <script_name>
#
#The slurm output file will by default, be written into the directory you submitted your job from  (slurm-JOBID.out)""".format(endpt,size))

    f.close()

    return


if __name__ == "__main__":
    main()
