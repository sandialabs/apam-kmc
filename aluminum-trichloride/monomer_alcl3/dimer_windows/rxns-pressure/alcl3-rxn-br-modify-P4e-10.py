#Quinn Campbell 1/21/2021  qcampbe@sandia.gov

#These rxns come from Radue et al 2021 https://doi.org/10.1021/acs.jpcc.1c00691


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

# -----------------------------------------------------------------------------
# AlCl3 adsorption
flux = 4e-10


elements_before.append(['u','u','u','*','*','*'])
elements_after.append(['AlCl2','u','Cl','*','*','*'])

kb  = 8.617e-5 # units of ev/K
rate_constant   =  FluxRate(flux)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# AlCl3 adsorption is not site specific as in the case of phosphine so let
# it adsorb on both sites

elements_before.append(['u','u','u','*','*','*'])
elements_after.append(['Cl','u','AlCl2','*','*','*'])

kb  = 8.617e-5 # units of ev/K
rate_constant   =  FluxRate(flux)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

elements_before.append(['Cl','u','Cl','*','*','*'])
elements_after.append(['u','u','u','*','*','*'])

kb  = 8.617e-5 # units of ev/K
rate_constant   = ArrheniusRate(4.85)

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))



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
#*********
#Now start the "double dimer" reactions
#*********
# # -----------------------------------------------------------------------------
# # AlCl3 adsorbs onto surface
#
# elements_before.append(['u','u','u','u','u','u'])
# elements_after.append(['AlCl2','u','AlCl2','Cl','u','Cl'])
#
# rate_constant   = FluxRate(flux)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))

# # -----------------------------------------------------------------------------
# # AlCl3 splits into 2
#
# elements_before.append(['AlCl2','u','AlCl2','u','u','u'])
# elements_after.append(['AlCl2','u','u','u','u','AlCl2'])
#
# rate_constant   = ArrheniusRate(1.26)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))
#
#
# # -----------------------------------------------------------------------------
# # The reverse reaction
#
# elements_before.append(['AlCl2','u','u','u','u','AlCl2'])
# elements_after.append(['AlCl2','u','AlCl2','u','u','u'])
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
#
# # -----------------------------------------------------------------------------
# # AlCl removes it's remaining Cl
#
# elements_before.append(['AlCl2','u','AlCl2','u','u','u'])
# elements_after.append(['AlCl','u','AlCl','Cl','u','Cl'])
#
# rate_constant   = ArrheniusRate(1.49)
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
# # the reverse
#
# elements_before.append(['AlCl','u','AlCl','Cl','u','Cl'])
# elements_after.append(['AlCl2','u','AlCl2','u','u','u'])
#
# rate_constant   = ArrheniusRate(0.22)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# AlCl2 loses Cl to nearby

elements_before.append(['AlCl2','u','*','u','u','*'])
elements_after.append(['AlCl','u','*','Cl','u','*'])

rate_constant   = ArrheniusRate(0.09)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# AlCl2 loses Cl to nearby

elements_before.append(['*','u','AlCl2','*','u','u'])
elements_after.append(['*','u','AlCl','*','u','Cl'])

rate_constant   = ArrheniusRate(0.09)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# AlCl2 straight to bridging AlCl

elements_before.append(['AlCl2','u','*','u','u','u'])
elements_after.append(['brAlCl','u','*','brAlCl','u','Cl'])

rate_constant   = ArrheniusRate(0.64)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# AlCl2 straight to bridging AlCl

elements_before.append(['*','u','AlCl2','u','u','u'])
elements_after.append(['*','u','brAlCl','Cl','u','brAlCl'])

rate_constant   = ArrheniusRate(0.64)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# AlCl2 straight to bridging AlCl

elements_before.append(['AlCl2','u','u','u','u','*'])
elements_after.append(['brAlCl','u','Cl','brAlCl','u','*'])

rate_constant   = ArrheniusRate(0.64)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# AlCl2 straight to bridging AlCl

elements_before.append(['u','u','AlCl2','*','u','u'])
elements_after.append(['Cl','u','brAlCl','*','u','brAlCl'])

rate_constant   = ArrheniusRate(0.64)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# AlCl to bridging AlCl

elements_before.append(['AlCl','u','*','Cl','u','u'])
elements_after.append(['brAlCl','u','*','brAlCl','u','Cl'])

rate_constant   = ArrheniusRate(0.56)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# AlCl to bridging AlCl

elements_before.append(['*','u','AlCl','u','u','Cl'])
elements_after.append(['*','u','brAlCl','Cl','u','brAlCl'])

rate_constant   = ArrheniusRate(0.56)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# # -----------------------------------------------------------------------------
# # bridging AlCl migration
#
# elements_before.append(['brAlCl','u','*','u','u','*'])
# elements_after.append(['u','u','*','brAlCl','u','*'])
#
# rate_constant   = ArrheniusRate(1.52)
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
# # bridging AlCl migration
#
# elements_before.append(['*','u','brAlCl','*','u','u'])
# elements_after.append(['*','u','u','*','u','brAlCl'])
#
# rate_constant   = ArrheniusRate(1.52)
# #print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# Cl migration

elements_before.append(['Cl','u','*','u','u','*'])
elements_after.append(['u','u','*','Cl','u','*'])

rate_constant   = ArrheniusRate(1.33)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# Cl migration

elements_before.append(['*','u','Cl','*','u','u'])
elements_after.append(['*','u','u','*','u','Cl'])

rate_constant   = ArrheniusRate(1.33)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# Cl diagonal migration

elements_before.append(['Cl','u','*','*','u','u'])
elements_after.append(['u','u','*','*','u','Cl'])

rate_constant   = ArrheniusRate(1.67)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# Cl diagonal migration

elements_before.append(['*','u','Cl','u','u','*'])
elements_after.append(['*','u','u','Cl','u','*'])

rate_constant   = ArrheniusRate(1.67)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))

# -----------------------------------------------------------------------------
# two alCl3 split up

elements_before.append(['AlCl3','u','AlCl3','u','u','u'])
elements_after.append(['AlCl2','u','Cl','Cl','u','AlCl2'])

rate_constant   = ArrheniusRate(0.43)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))


# -----------------------------------------------------------------------------
# two alCl3 split up

elements_before.append(['AlCl3','u','AlCl3','u','u','u'])
elements_after.append(['Cl','u','AlCl2','AlCl2','u','Cl'])

rate_constant   = ArrheniusRate(0.43)
#print rate_constant

processes.append(KMCProcess(
    coordinates=coordinates,
    elements_before=elements_before[-1],
    elements_after=elements_after[-1],
    basis_sites=basis_sites,
    rate_constant=rate_constant))


# # -----------------------------------------------------------------------------
# # hydrogen hoppping
#
# elements_before.append(['H','*','*','u','*','*'])
# elements_after.append(['u','*','*','H','*','*'])
#
# #this barrier subject to change!
# rate_constant   = ArrheniusRate(1.51)
# print rate_constant
#
# processes.append(KMCProcess(
#     coordinates=coordinates,
#     elements_before=elements_before[-1],
#     elements_after=elements_after[-1],
#     basis_sites=basis_sites,
#     rate_constant=rate_constant))



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


interactions = KMCInteractions(
    processes=processes,
    implicit_wildcards=True)

