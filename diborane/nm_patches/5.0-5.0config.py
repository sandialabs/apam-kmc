#
# This is the lattice setup for a silicon surface

from KMCLib import *
import math

# ----------------------------------------------------------------------------
# Unit cell

cell_vectors = [[   0.3874e+00*2.,   0.000000e+00,   0.000000e+00],
                [   0.000000e+00,   0.7666e+00,   0.000000e+00],
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

#in units of nm
height = 5.0
width = 5.0

dimer_height = 0.766
dimer_width = 0.3874

num_dimers_wide = int(math.ceil(width/dimer_width))
num_dimers_tall = int(math.ceil(height/dimer_height))

reps_x = int(math.ceil((num_dimers_wide+3.)/2.))
reps_y = int(math.ceil((num_dimers_tall+3.)))
#print reps

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
    if (site[0] < 1.0) or (site[0] >= (num_dimers_wide +2)*0.5) or (site[1]<1.0) or (site[1] >= (num_dimers_tall+2.)):
        ending_decimal = site[1]-math.floor(site[1])
        #print ending_decimal
        if abs(ending_decimal - 0.333333) < 0.1 or abs(ending_decimal - 1.333333) < 0.1:
            types.append('u')
        else:
            types.append('H')
    else:
        types.append('u')

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
