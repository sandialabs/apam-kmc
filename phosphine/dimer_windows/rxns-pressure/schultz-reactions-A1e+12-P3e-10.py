#Quinn Campbell 6/21/2019  qcampbe@sandia.gov

#These rxns come from Warschkow 2016 and Peter Schultz 2017 MRS talk

from KMCLib import *


def ArrheniusRate(barrier):
    #For now setting these here, Ideally we'd have a way to pass in this variable.
    A   = 1.0e12
    T   = 300
    kb  = 8.617e-5 # units of ev/K
    from numpy import exp as np_exp
    return A*np_exp(-barrier/(kb*T))



def FluxRate(pressure,M=33.99758,T=300):
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

#Adsorption of PH3
#We will continue the tradition of assuming roughly 1/s

coordinates = [[   0.000000e+00,   0.000000e+00,   0.000000e+00],
           [   0.000000e+00,   0.333333e+00,   0.000000e+00],
           [   0.000000e+00,   0.666666e+00,   0.000000e+00],
           [   0.500000e+00,   0.000000e+00,   0.000000e+00],
           [   0.500000e+00,   0.333333e+00,   0.000000e+00],
           [   0.500000e+00,   0.666666e+00,   0.000000e+00]]

basis_sites     = [0]

#Am going to be fairly liberal on my use of wildcards for now. If it comes back
# to haunt me, then you heard it here first folks, wildcards are the devil

elements_before.append(['u','u','*','*','*','*'])
elements_after.append(['PH3','u','*','*','*','*'])

#rate_constant   = FluxRate(1e-6)
#The Fuschle data uses the pressure of 3.75e-8
rate_constant   = FluxRate(3e-10)

print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# Desorption energy
# Rxn 0

elements_before.append(['PH3','u','*','*','*','*'])
elements_after.append(['u','u','*','*','*','*'])

rate_constant   = ArrheniusRate(0.91)
#print rate_constant

print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# lose H to same dimer row
# Rxn 1

elements_before.append(['PH3','u','u','*','*','*'])
elements_after.append(['PH2','u','H','*','*','*'])

rate_constant   = ArrheniusRate(0.62)
#print rate_constant

print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# Warschkow B2 --> B4
# Rxn 2

elements_before.append(['PH2','u','u','*','*','*'])
elements_after.append(['b','PH2','b','*','*','*'])

rate_constant   = ArrheniusRate(0.76)
#print 'rate constant: ' , rate_constant
print 1.0/rate_constant


processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# Wilson fig. 3c

elements_before.append(['PH2','u','H','*','*','*'])
elements_after.append(['u','u','u','*','*','*'])

rate_constant   = ArrheniusRate(2.0)
#print rate_constant
print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))

# # -----------------------------------------------------------------------------
# # Wilson fig. 3c
#
# elements_before.append(['H','u','PH2','*','*','*'])
# elements_after.append(['u','u','u','*','*','*'])
#
# rate_constant   = ArrheniusRate(2.0)
# #print rate_constant
# print 1.0/rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# Warschkow B4 ---> B2

elements_before.append(['b','PH2','b','*','*','*'])
elements_after.append(['PH2','u','u','*','*','*'])

rate_constant   = ArrheniusRate(0.13)
print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# Warschkow B4 ---> B3

elements_before.append(['b','PH2','b','*','*','*'])
elements_after.append(['u','u','PH2','*','*','*'])

rate_constant   = ArrheniusRate(0.09)
print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# Warschkow B3 ---> B4

elements_before.append(['u','u','PH2','*','*','*'])
elements_after.append(['b','PH2','b','*','*','*'])

rate_constant   = ArrheniusRate(0.73)
print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# # Wilson fig. 3c
#
# elements_before.append(['H','u','H','*','*','*'])
# elements_after.append(['u','u','u','*','*','*'])
#
# rate_constant   = ArrheniusRate(2.1)
# #print rate_constant
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
# -----------------------------------------------------------------------------
# H migration to next dimer
# Rxn 3

elements_before.append(['PH3','u','*','u','u','*'])
elements_after.append(['PH2','u','*','H','u','*'])

rate_constant   = ArrheniusRate(0.35)
#print rate_constant

print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# Warschkow B4 ---> C1
# Rxn 4

elements_before.append(['b','PH2','b','u','u','*'])
elements_after.append(['b','PH','b','H','u','*'])

rate_constant   = ArrheniusRate(0.36)
#print rate_constant

print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# Warschkow B4 ---> C1 the alternative throws off the up down pairing but
# should be fine.
# Rxn 4

elements_before.append(['b','PH2','b','*','u','u'])
elements_after.append(['b','PH','b','*','u','H'])

rate_constant   = ArrheniusRate(0.36)
#print rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))





# -----------------------------------------------------------------------------
# Warschkow B1 ---> B3
# Rxn 5

elements_before.append(['PH2','u','*','u','u','*'])
elements_after.append(['u','u','*','PH2','u','*'])

rate_constant   = ArrheniusRate(0.92)
#print rate_constant
print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# Warschkow B3 ---> B1
# Rxn 5

elements_before.append(['u','u','*','PH2','u','*'])
elements_after.append(['PH2','u','*','u','u','*'])

rate_constant   = ArrheniusRate(0.68)
#print rate_constant
print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# Warschkow B3 ---> C1
# Rxn 6

elements_before.append(['u','u','*','PH2','u','u'])
elements_after.append(['H','u','*','b','PH','b'])

rate_constant   = ArrheniusRate(0.74)
#print 'rate constant: ' , rate_constant
print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# Warschkow B3 ---> C1 (alternate H version, since H don't care about your up/down
# PH3 ordering)
# Rxn 6

elements_before.append(['*','u','u','PH2','u','u'])
elements_after.append(['*','u','H','b','PH','b'])

rate_constant   = ArrheniusRate(0.74)
#print rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))

# # -----------------------------------------------------------------------------
# # Warschkow B3 ---> C1
# # Rxn 6
#
# elements_before.append(['u','u','*','u','u','PH2'])
# elements_after.append(['H','u','*','b','PH','b'])
#
# rate_constant   = ArrheniusRate(0.74)
# #print 'rate constant: ' , rate_constant
# #print 1.0/rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))
#
# # -----------------------------------------------------------------------------
# # Warschkow B3 ---> C1 (alternate H version, since H don't care about your up/down
# # PH3 ordering)
# # Rxn 6
#
# elements_before.append(['*','u','u','u','u','PH2'])
# elements_after.append(['*','u','H','b','PH','b'])
#
# rate_constant   = ArrheniusRate(0.74)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# Concerted reaction
# Rxn 7

elements_before.append(['PH2','u','u','u','*','*'])
elements_after.append(['b','PH','b','H','*','*'])

rate_constant   = ArrheniusRate(0.96)
#print rate_constant
print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# Wilson fig 4c
# Rxn 7

elements_before.append(['H','u','u','H','u','u'])
elements_after.append(['u','u','u','u','u','u'])

rate_constant   = ArrheniusRate(1.7)
#print rate_constant
print 1.0/rate_constant

processes.append(KMCProcess(
coordinates=coordinates,
elements_before=elements_before[-1],
elements_after=elements_after[-1],
basis_sites=basis_sites,
rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# Wilson fig 3e
# Rxn 7

elements_before.append(['PH2','u','*','H','u','*'])
elements_after.append(['u','u','*','u','u','*'])

rate_constant   = ArrheniusRate(2.0)
#print rate_constant
print 1.0/rate_constant

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

# for i in range(len(processes)):
#     print i
#     print processes[i]._coordinates
#     print processes[i].elementsBefore()
#     rough_time = 1.0 / processes[i].rateConstant()
#     print rough_time
#     #print processes[i].rateConstant()
#     print processes[i].elementsAfter()
#
# print len(processes)


interactions = KMCInteractions(
processes=processes,
implicit_wildcards=True)
