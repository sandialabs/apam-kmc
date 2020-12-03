#Quinn Campbell  qcampbe@sandia.gov

from KMCLib import *


def ArrheniusRate(barrier):
    #For now setting these here, Ideally we'd have a way to pass in this variable.
    A   = 1.0e12
    T   = 393 # default is 393
    kb  = 8.617e-5 # units of ev/K
    from numpy import exp as np_exp
    return A*np_exp(-barrier/(kb*T))



def FluxRate(pressure,M=27.668,T=393):
    dimerArea_in_cm2 = 0.52918**2*1e-14
    torr2pa = 133.322368
    from numpy import sqrt as np_sqrt
    return (2.63e20*pressure*torr2pa/np_sqrt(M*T))*dimerArea_in_cm2


#since the coordinate order gets messed up, if you want to do a straightforward
# reversal, it makes more sense to do it on your coordinates than on the ones
# you could retrieve afterwards as I was doing. Somewhat cumbersome, but much
# more straightforward for the coder (ahem, yours truly) to wrap their mind around
processes = []
elements_before = []
elements_after = []

# -----------------------------------------------------------------------------
# Interactions

#As in previous code, it makes most sense to add all reactions on a single dimer
#at once and then to swtich to "double dimer reactions" so that the symmetrry
# can be maintained


coordinates = [[   0.000000e+00,   0.000000e+00,   0.000000e+00],
               [   0.000000e+00,   0.333333e+00,   0.000000e+00],
               [   0.000000e+00,   0.666666e+00,   0.000000e+00],
               [   0.500000e+00,   0.000000e+00,   0.000000e+00],
               [   0.500000e+00,   0.333333e+00,   0.000000e+00],
               [   0.500000e+00,   0.666666e+00,   0.000000e+00]]

basis_sites     = [0]


# -----------------------------------------------------------------------------
# diborane adsorption

elements_before.append(['u','u','u','*','*','*'])
elements_after.append(['B2H6','u','u','*','*','*'])

T   = 300
kb  = 8.617e-5 # units of ev/K
barrier = 0.32
from numpy import exp as np_exp
rate_constant   =FluxRate(1.5e-7)#*np_exp(-barrier/(kb*T))
print "this is the rate constant ", rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# diborane splitting

elements_before.append(['B2H6','u','u','*','*','*'])
elements_after.append(['BH3','u','BH3','*','*','*'])


barrier = 0.32
from numpy import exp as np_exp
rate_constant = ArrheniusRate(barrier)
print "this is the rate constant ", rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# diborane desorption

elements_before.append(['BH3','u','BH3','*','*','*'])
elements_after.append(['u','u','u','*','*','*'])

rate_constant   =ArrheniusRate(1.41)

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH3 desorption

elements_before.append(['BH3','u','u','*','*','*'])
elements_after.append(['u','u','u','*','*','*'])

rate_constant   =ArrheniusRate(1.21)

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH3 desorption

elements_before.append(['u','u','BH3','*','*','*'])
elements_after.append(['u','u','u','*','*','*'])

rate_constant   =ArrheniusRate(1.21)

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH3 desorption

elements_before.append(['BH3','u','H','*','*','*'])
elements_after.append(['u','u','H','*','*','*'])

rate_constant   =ArrheniusRate(1.21)

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# BH3 desorption

elements_before.append(['H','u','BH3','*','*','*'])
elements_after.append(['H','u','u','*','*','*'])

rate_constant   =ArrheniusRate(1.21)

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))



# -----------------------------------------------------------------------------
# br BH2 + H desorption

elements_before.append(['u','BH3','u','*','*','*'])
elements_after.append(['u','u','u','*','*','*'])

rate_constant   =ArrheniusRate(1.9)

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# # -----------------------------------------------------------------------------
# # br BH2 + H desorption
#
# elements_before.append(['*','BH2','H','*','*','*'])
# elements_after.append(['*','u','u','*','*','*'])
#
# rate_constant   =ArrheniusRate(1.9)
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# H2 desorption

elements_before.append(['H','u','H','*','*','*'])
elements_after.append(['u','u','u','*','*','*'])

rate_constant   =ArrheniusRate(2.1)

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))


# # diborane adsorption
#
# elements_before.append(['*','*','*','u','u','u'])
# elements_after.append(['*','*','*','BH3','u','BH3'])
#
# rate_constant   = ArrheniusRate(0.32)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH3 moving to middle position

elements_before.append(['BH3','u','u','*','*','*'])
elements_after.append(['u','BH3','u','*','*','*'])

rate_constant   = ArrheniusRate(0.04)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH2 moving to middle position

elements_before.append(['BH2','u','u','*','*','*'])
elements_after.append(['u','BH2','u','*','*','*'])

rate_constant   = ArrheniusRate(0.008)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH3 moving to middle position

elements_before.append(['u','u','BH3','*','*','*'])
elements_after.append(['u','BH3','u','*','*','*'])

rate_constant   = ArrheniusRate(0.04)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH2 moving to middle position

elements_before.append(['u','u','BH2','*','*','*'])
elements_after.append(['u','BH2','u','*','*','*'])

rate_constant   = ArrheniusRate(0.008)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))


# # -----------------------------------------------------------------------------
# # alternate two H desorption
#
# elements_before.append(['H','u','H','*','*','*'])
# elements_after.append(['u','u','u','*','*','*'])
#
# #rate_constant   = ArrheniusRate(1.69)
# rate_constant   = ArrheniusRate(2.1)
# #print rate_constant
#
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))


sym_up_to = len(processes)
#Only one direction of symmetry here to avoid double counting of possibilities
# this is the distinction between the "single" and "double" dimer processes in the
# previous codes. Both are on a double dimer lattice now, but the "single" dimer
# reactions only get checked in one direction to preserve up down sym, whereas
# double dimer reactions get checked in both directions
for i in range(sym_up_to):
    elements_b = elements_before[i]
    elements_a = elements_after[i]
    #print elements_a
    elements_b =elements_b[3:][::-1]+ elements_b[:3][::-1]
    elements_a = elements_a[3:][::-1]+ elements_a[:3][::-1]
    #print elements_a
    elements_before.append(elements_b)
    elements_after.append(elements_a)
    #print elements_before[-1]
    #print elements_after[-1]
    processes.append(KMCProcess(
        coordinates=coordinates,
        elements_before=elements_before[-1],
        elements_after=elements_after[-1],
        basis_sites=basis_sites,
        rate_constant=processes[i].rateConstant()
        ))

sym_start = len(processes)

#Now start the "double dimer" reactions

# -----------------------------------------------------------------------------
# diborane splitting into 2 BH3 (the good outcome)

elements_before.append(['BH3','u','BH3','u','u','u'])
elements_after.append(['u','BH3','u','u','u','BH3'])

rate_constant   = ArrheniusRate(0.91)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# diborane splitting into 2 BH3 (the good outcome)

elements_before.append(['BH3','u','BH3','u','u','u'])
elements_after.append(['u','BH3','u','BH3','u','u'])

rate_constant   = ArrheniusRate(0.91)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# reverse rxn: diborane splitting into 2 BH3 (the good outcome)
#

elements_before.append(['u','BH3','u','u','u','BH3'])
elements_after.append(['BH3','u','BH3','u','u','u'])

rate_constant   = ArrheniusRate(0.99)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# reverse rxn: diborane splitting into 2 BH3 (the good outcome)
#

elements_before.append(['u','BH3','u','BH3','u','u'])
elements_after.append(['BH3','u','BH3','u','u','u'])

rate_constant   = ArrheniusRate(0.99)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# # -----------------------------------------------------------------------------
# # diborane losing a hydrogen. The *slightly* more likely outcome
#
# elements_before.append(['BH3','u','BH3','u','u','u'])
# elements_after.append(['BH3','u','BH2','u','u','H'])
#
# rate_constant   = ArrheniusRate(0.89)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))
#
# # -----------------------------------------------------------------------------
# # reverse rxn: diborane losing a hydrogen. The *slightly* more likely outcome
#
# elements_before.append(['BH3','u','BH2','u','u','H'])
# elements_after.append(['BH3','u','BH3','u','u','u'])
#
# rate_constant   = ArrheniusRate(1.57)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# diborane losing a hydrogen. The *slightly* more likely outcome

elements_before.append(['BH3','u','BH3','u','u','u'])
elements_after.append(['BH2','u','BH3','H','u','u'])

rate_constant   = ArrheniusRate(0.89)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# reverse rxn: diborane losing a hydrogen. The *slightly* more likely outcome

elements_before.append(['BH2','u','BH3','H','u','u'])
elements_after.append(['BH3','u','BH3','u','u','u'])

rate_constant   = ArrheniusRate(1.57)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))


# # -----------------------------------------------------------------------------
# # losing the second hydrogen. it's gonna happen if you're here
#
# elements_before.append(['BH3','u','BH2','u','u','H'])
# elements_after.append(['BH2','u','BH2','H','u','H'])
#
# rate_constant   = ArrheniusRate(0.11)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))


# # -----------------------------------------------------------------------------
# # reverse: losing the second hydrogen. it's gonna happen if you're here
#
# elements_before.append(['BH2','u','BH2','H','u','H'])
# elements_after.append(['BH3','u','BH2','u','u','H'])
#
# rate_constant   = ArrheniusRate(1.34)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# losing the second hydrogen. it's gonna happen if you're here

elements_before.append(['BH2','u','BH3','H','u','u'])
elements_after.append(['BH2','u','BH2','H','u','H'])

rate_constant   = ArrheniusRate(0.11)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))


# # -----------------------------------------------------------------------------
# # reverse: losing the second hydrogen. it's gonna happen if you're here
#
# elements_before.append(['BH2','u','BH2','H','u','H'])
# elements_after.append(['BH2','u','BH3','H','u','u'])
#
# rate_constant   = ArrheniusRate(1.34)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# diborane losing a hydrogen. The *slightly* more likely outcome

elements_before.append(['BH2','u','BH2','u','u','u'])
elements_after.append(['u','BH2','u','u','BH2','u'])

rate_constant   = ArrheniusRate(1.25)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# reverse: BH2 splitting

elements_before.append(['u','BH2','u','u','BH2','u'])
elements_after.append(['BH2','u','BH2','u','u','u'])

rate_constant   = ArrheniusRate(1.33)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# losing H to nearby dimer

elements_before.append(['BH2','u','BH2','u','u','u'])
elements_after.append(['BH','u','BH2','H','u','u'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.22)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# reverse: losing H to nearby dimer

elements_before.append(['BH','u','BH2','H','u','u'])
elements_after.append(['BH2','u','BH2','u','u','u'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(0.67)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# # -----------------------------------------------------------------------------
# # losing H to nearby dimer
#
# elements_before.append(['BH2','u','BH2','u','u','u'])
# elements_after.append(['BH2','u','BH','u','u','H'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(1.22)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))
#
# # -----------------------------------------------------------------------------
# # reverse: losing H to nearby dimer
#
# elements_before.append(['BH2','u','BH','u','u','H'])
# elements_after.append(['BH2','u','BH2','u','u','u'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(0.67)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))





# -----------------------------------------------------------------------------
# losing the second hydrogen

elements_before.append(['BH','u','BH2','H','u','u'])
elements_after.append(['BH','u','BH','H','u','H'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(0.01)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))


# # -----------------------------------------------------------------------------
# # reverse: losing the second hydrogen
#
# elements_before.append(['BH','u','BH','H','u','H'])
# elements_after.append(['BH','u','BH2','H','u','u'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(1.10)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))

# # -----------------------------------------------------------------------------
# # losing the second hydrogen
#
# elements_before.append(['BH2','u','BH','u','u','H'])
# elements_after.append(['BH','u','BH','H','u','H'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(0.01)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))


# # -----------------------------------------------------------------------------
# # reverse: losing the second hydrogen
#
# elements_before.append(['BH','u','BH','H','u','H'])
# elements_after.append(['BH','u','BH2','H','u','u'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(1.10)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))




# -----------------------------------------------------------------------------
# BH3 moving to  bridging position, losing H to nearby row

elements_before.append(['BH3','u','u','u','u','*'])
elements_after.append(['u','BH2','u','H','u','*'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(0.82)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# reverse: BH3 moving to  bridging position, losing H to nearby row

elements_before.append(['u','BH2','u','H','u','*'])
elements_after.append(['BH3','u','u','u','u','*'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.84)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# hydrogen hoppping

elements_before.append(['H','u','u','u','u','u'])
elements_after.append(['u','u','u','H','u','u'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.73)
print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# hydrogen hoppping

elements_before.append(['u','u','H','u','u','u'])
elements_after.append(['u','u','u','u','u','H'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.73)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# BH2 losing a hydrogen

elements_before.append(['*','BH2','*','u','u','H'])
elements_after.append(['*','brBH','*','H','u','H'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.34)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH2 losing a hydrogen

elements_before.append(['*','BH2','*','H','u','u'])
elements_after.append(['*','brBH','*','H','u','H'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.34)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))
#
# # -----------------------------------------------------------------------------
# # BH2 losing a hydrogen
#
# elements_before.append(['*','brBH','*','H','u','*'])
# elements_after.append(['*','BH2','*','u','u','*'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(1.41)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))
#
# # -----------------------------------------------------------------------------
# # BH2 losing a hydrogen
#
# elements_before.append(['*','brBH','*','*','u','H'])
# elements_after.append(['*','BH2','*','*','u','u'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(1.41)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH3 losing a hydrogen

elements_before.append(['*','BH3','*','u','u','u'])
elements_after.append(['*','BH2','*','H','u','u'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.51)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH3 losing a hydrogen

elements_before.append(['*','BH3','*','u','u','u'])
elements_after.append(['*','BH2','*','u','u','H'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.51)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH3 losing a hydrogen

elements_before.append(['*','BH2','*','H','u','u'])
elements_after.append(['*','BH3','*','u','u','u'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.83)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH3 losing a hydrogen

elements_before.append(['*','BH2','*','u','u','H'])
elements_after.append(['*','BH3','*','u','u','u'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.83)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))


# # -----------------------------------------------------------------------------
# # BH losing a hydrogen
#
# elements_before.append(['*','BH','*','u','*','*'])
# elements_after.append(['*','B','*','H','*','*'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(1.89)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))
#
# # -----------------------------------------------------------------------------
# # BH losing a hydrogen
#
# elements_before.append(['*','BH','*','*','*','u'])
# elements_after.append(['*','B','*','*','*','H'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(1.89)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# BH2 hoppping dimers

elements_before.append(['u','BH2','u','u','u','u'])
elements_after.append(['u','u','u','u','BH2','u'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.63)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# BH2 hoppping dimers

elements_before.append(['u','u','u','u','BH2','u'])
elements_after.append(['u','BH2','u','u','u','u'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.63)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# # -----------------------------------------------------------------------------
# # BH moving to a bridging position between two dimers
#
# elements_before.append(['u','BH','u','*','*','u'])
# elements_after.append(['u','u','brBH','*','*','brBH'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(1.29)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))
#
# # -----------------------------------------------------------------------------
# # BH moving to a bridging position between two dimers
#
# elements_before.append(['u','BH','u','u','*','*'])
# elements_after.append(['brBH','u','u','brBH','*','*'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(1.29)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))
#
# # -----------------------------------------------------------------------------
# # BH moving to a bridging position between two dimers
#
# elements_before.append(['u','u','brBH','*','*','brBH'])
# elements_after.append(['u','BH','u','*','*','u'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(1.55)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))
#
# # -----------------------------------------------------------------------------
# # BH moving to a bridging position between two dimers
#
# elements_before.append(['brBH','u','u','brBH','*','*'])
# elements_after.append(['u','BH','u','u','*','*'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(1.55)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# BH moving to a bridging position between two dimers

elements_before.append(['H','u','u','H','u','u'])
elements_after.append(['u','u','u','u','u','u'])

#this barrier subject to change!
rate_constant   = ArrheniusRate(1.7)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

sym_up_to = len(processes)

for i in range(sym_start,sym_up_to):
    elements_b = elements_before[i]
    elements_a = elements_after[i]
    elements_b =elements_b[3:][::-1]+ elements_b[:3][::-1]
    elements_a = elements_a[3:][::-1]+ elements_a[:3][::-1]
    elements_before.append(elements_b)
    elements_after.append(elements_a)
    processes.append(KMCProcess(
        coordinates=coordinates,
        elements_before=elements_before[-1],
        elements_after=elements_after[-1],
        basis_sites=basis_sites,
        rate_constant=processes[i].rateConstant()
        ))
    #print elements_b

sym_up_to = len(processes)

new_coords = [[   0.000000e+00,   0.000000e+00,   0.000000e+00],
               [   0.000000e+00,   0.333333e+00,   0.000000e+00],
               [   0.000000e+00,   0.666666e+00,   0.000000e+00],
               [   -0.500000e+00,   0.000000e+00,   0.000000e+00],
               [   -0.500000e+00,   0.333333e+00,   0.000000e+00],
               [   -0.500000e+00,   0.666666e+00,   0.000000e+00]]

for i in range(sym_start,sym_up_to):
    elements_before.append(elements_before[i])
    elements_after.append(elements_after[i])
    processes.append(KMCProcess(
        coordinates=new_coords,
        elements_before=elements_before[-1],
        elements_after=elements_after[-1],
        basis_sites=basis_sites,
        rate_constant=processes[i].rateConstant()
        ))




#NOTE: the internal KMCprocess method reorders the coordinates so
# the output will look wrong. Thus check coordinate order before freaking out
"""
for i in range(len(processes)):
    print processes[i]._coordinates
    print processes[i].elementsBefore()
    print processes[i].rateConstant()
    print processes[i].elementsAfter()

print len(processes)
"""

interactions = KMCInteractions(
    processes=processes,
    implicit_wildcards=True)
