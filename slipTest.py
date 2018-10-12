# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 16:40:38 2018

@author: bmusammartanoV
"""
import numpy as num
import matplotlib.pyplot as pl

DPi = 700 #Pressure drop
Qi= 0.5  #Flow rate cm/s
rot = 1.2 #density
slip = 0.64 #slip factor

#Impeller Inlet Geometry
B1 = 0.05  #inlet width
D1 = 0.207 #inner diam
w = 157.07 #ras/s
Nb = 10 #Number of blades
th = 0.002 #blades thickness

################
##Slip Formulas
###############

###############
#Stanitz's
i = []
S1 = []
c = 0
corr = 999
toll1 = 1e-20
while corr > toll1:
    Pow = 1.1*DPi*Qi   #Effective power
    Ws  = Pow/(rot*Qi) #Work to be done
    U2  = num.sqrt(1.1*DPi/(rot*slip)) #Power = Peuler = rho*Q*(U2*0.8*U2)
    D2  = 2 * U2/w #impeller outer diameter
    A2 = 0.25*num.pi*D2**2 #impeller outer aerea
    B2 = B1

    V2m = Qi/((num.pi*D2 - Nb*th)*B2) #
    V2t = slip*U2
    V2 = num.sqrt(V2m**2 + V2t**2)
    alpha2 = num.rad2deg(num.arctan(V2m/V2t)) #alpha2 in deg

    W2t = (1-slip) * U2
    W2 = num.sqrt( V2m**2 + W2t**2 )
    beta2 = num.arctan(V2m/W2t) #beta2 in deg
    phi2 = (W2/U2)
    if beta2 < 45:
        beta2 = 45
    #Stanitz's
    slipN = 1 - (1.98/Nb)/(1-phi2/num.tan(beta2))
    corr = abs(slipN - slip)
    i.append(c)
    c += 1
    slip = slipN
    S1.append(slip)

    print("Beta2 = {:3.4f} deg".format((beta2)))
    print("{}th slip factor estimation: {:3.4f} ".format(c,slip))
################
    
################
#Balje's formula
i2 =[]
S2 = []
c = 0
toll2 = 1e-20
corr=999
while corr > toll2:
    Pow = 1.1*DPi*Qi   #Effective power
    Ws  = Pow/(rot*Qi) #Work to be done
    U2  = num.sqrt(1.1*DPi/(rot*slip)) #Power = Peuler = rho*Q*(U2*0.8*U2)
    D2  = 2 * U2/w #impeller outer diameter
    A2 = 0.25*num.pi*D2**2 #impeller outer aerea
    B2 = B1

    V2m = Qi/((num.pi*D2 - Nb*th)*B2) #
    V2t = slip*U2
    V2 = num.sqrt(V2m**2 + V2t**2)
    alpha2 = num.rad2deg(num.arctan(V2m/V2t)) #alpha2 in deg

    W2t = (1-slip) * U2
    W2 = num.sqrt( V2m**2 + W2t**2 )
    beta2 = num.arctan(V2m/W2t) #beta2 in deg
    phiD = D1/D2
    #Balje's formula
    slipN = (1+6.2/(Nb*phiD**(2/3)))**-1
    corr = abs(slipN - slip)
    i2.append(c)
    c += 1
    slip = slipN
    S2.append(slip)

    print("Beta2 = {:3.4f} deg".format(num.rad2deg(beta2)))
    print("{}th slip factor estimation: {:3.4f} ".format(c,slip))

pl.plot(i,S1,"-o",label="Stanitz's formula")
pl.plot(i2,S2,"-*",label="Balje's formula")
pl.legend()

pl.title("Slip Factor - step {}".format(c))
pl.show()


#Library
#https://en.wikipedia.org/wiki/Slip_factor