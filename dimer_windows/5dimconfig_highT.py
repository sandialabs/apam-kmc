#
# This is the lattice setup for a silicon surface

from KMCLib import *
import math

# ----------------------------------------------------------------------------
# Unit cell

cell_vectors = [[   2.000000e+00,   0.000000e+00,   0.000000e+00],
                [   0.000000e+00,   3.000000e+00,   0.000000e+00],
                [   0.000000e+00,   0.000000e+00,   1.000000e+00]]

#Here 0.00 refers to the first dimer site
# y=1 is the bridging site (if a PH gets here, we will treat it as incorp, otherwise mostly ignored)
# y = 2 is the second dimer site
basis_points = [[   0.000000e+00,   0.000000e+00,   0.000000e+00],
                [   0.000000e+00,   0.333333e+00,   0.000000e+00],
                [   0.000000e+00,   0.666666e+00,   0.000000e+00],
                [   0.500000e+00,   0.000000e+00,   0.000000e+00],
                [   0.500000e+00,   0.333333e+00,   0.000000e+00],
                [   0.500000e+00,   0.666666e+00,   0.000000e+00],]



unit_cell = KMCUnitCell(
    cell_vectors=cell_vectors,
    basis_points=basis_points)

# -----------------------------------------------------------------------------
# Lattice

num_dimers = 5

reps = int(math.ceil((num_dimers+3.)/2.))
#print reps

lattice = KMCLattice(
    unit_cell=unit_cell,
    repetitions=(reps,1,1),
    periodic=(False, False, False))

#Starting with a totally unoccupied cell
trajfile="5dimwindow.py"
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
        print final_index
        break


types = elem[final_index]
#print types

#I will use 'u' for unoccupied
# and 'b' for involved in a bridging reaction
# otherwise types will match their atomic label (i.e PH2)
#I'm including a few types such as H2 and P that we essentially won't run into
#Shouldn't be an issue, but if it's causing issues later, here's your problem
possible_types = ['u','brBH','BH3','BH2','BH','B','H','H2','B2H6']

configuration = KMCConfiguration(
    lattice=lattice,
    types=types,
    possible_types=possible_types)
