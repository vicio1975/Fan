# -*- coding: utf-8 -*-
"""
Created on Wed March 14 09:21:55 2018
Version 1

Goals:
        1. Design Fan geometry - incompressible flow
            1.1 Write a file with the design parameters
        2. Design Fan geometry - compressible flow
            2.1 Write a file with the design parameters
        3. Different blades configurations (radial tip, sheet and aerofoil blades ...)
        4. Automatic CAD

Note1 >> For Boundary layer use the tool BoundaryLayer.py !!

@author: Vincenzo Sammartano

"""

### Libraries
import numpy as num
from math import pi
#import matplotlib.pyplot as plt
#import sys
#import shutil

### Functions

def impellerIn(Qi):
    #inlet Duct and impeller inlet design
    #Duct and impeller parameters
    teta = num.radians(30) # duct contraction
    ind = 0.2  #induction improvement
    uvRatio = 1.1 #ratio U/V
    #Duct geometry and Impeller eye geometry
    D1 = num.power((uvRatio * 8/pi) * (Qi/w),1/3)
    Dduct = (1+ind)*D1

    A1 = 0.25 * pi * D1**2
    Aduct = 0.25 * pi * Dduct**2
    Lduct = 0.1 * D1 / (2 * num.tan(teta))
    #Velocity triangle
    U1 = w * D1/2
    V1 = U1/1.1
    
    Vduct = 4*Qi/(pi*(Dduct**2))
    W1 = num.power(U1**2+V1**2,0.5)
    beta1 = num.rad2deg(num.arctan(1/1.1))
    #impeller width at inlet section 
    B1 = (Qi/V1)/(pi*D1 - Nb*th)
    return [D1,Dduct,Lduct,U1,V1,Vduct,W1,beta1,B1,A1,Aduct]

def impellerOut(DPi,Qi,rot,slip,B1): #qL,DPi,Qi,rot,slip,B1,beta2
    Pow = 1.1*abs(DPi)*Qi   #Effective power
    Ws  = Pow/(rot*Qi) #Work to be done
    U2  = num.sqrt(abs(DPi) * 1.1/(rot*slip)) #Power = Peuler = rho*Q*(U2*0.8*U2)
    D2  = 2 * U2/w #impeller outer diameter
    A2 = 0.25*pi*D2**2 #impeller outer aerea
    B2 = B1

    V2m = Qi/((pi*D2 - Nb*th)*B2) #meridian velocity component
    V2t = slip*U2 #tangential velocity component
    V2 = num.sqrt(V2m**2 + V2t**2) #outlet velocity800
    alpha2 = num.rad2deg(num.arctan(V2m/V2t)) #alpha2 in deg

    WU2 = (1-slip) * U2
    W2 = num.sqrt( V2m**2 + WU2**2 )
    beta2 = num.rad2deg(num.arctan(V2m/WU2)) #beta2 in deg
    return [U2,D2,V2,V2m,W2,alpha2,beta2,B2,A2,Pow,Ws,slip]

def blades(D1,D2,beta1):
    #modify in case of beta2<90
    Rb = 0.5*((D2/2)**2-(D1/2)**2)/((D1/2)*num.cos(num.radians(beta1)))
    return Rb

def volute(V1,Ws,rot,DPi,Qi,B1,D2,alpha2):
    clv = 0.015 * D2 # clearance = 1.5% D2
    V4 = num.sqrt(V1**2 + 2*Ws - 2*abs(DPi)/rot)
    D3 = D2 + 2*clv
    R3 = D3/2
    Bv = 2.5 * B1
    R4 = Qi/(V4*Bv) + R3
    DR = R4 - R3
    R2 = D2/2 #Radius of the outer impeller
    Rt = R3   #Radius of the volute tongue
    tetaT = 132*num.log10(Rt/R2)/num.tan(num.radians(alpha2))
    startAng = tetaT
    stopAng = 360
    teta = num.arange(startAng,stopAng,2)
    Rtet = Rt + (teta/360) * DR
    Xtet = Rtet * num.sin(num.radians(teta))
    Ytet = Rtet * num.cos(num.radians(teta))
    return [V4,D3,Bv,R4,DR,Xtet,Ytet]

def design(nameF,Qi,Pin,DPi,Nb,th,w,T):
    Rf = 287.058 # Universal Constant of Gases [J/(Kg K)]
    pAtm = 101325 # [Pa] atmospheric pressure
    slip = 1-(1.98/Nb) #Stanitz
    #Taking compressibility effect into consideration
    rot = pAtm/(Rf*T) #air density at T celsius
    corr = False
    QL = 0 #Flow rate

    while corr == False:
        #design with the input data Q & DP
        D1,Dduct,Lduct,U1,V1,Vduct,W1,beta1,B1,A1,Aduct = impellerIn(Qi)
        #function to estimate the outlet impeller geometry  - F(Qi,DPi)
        U2,D2,V2,V2m,W2,alpha2,beta2,B2,A2,Pow,Ws,slipf = impellerOut(DPi,Qi,rot,slip,B1)
        #function to estimate the volute casing geometry - F(Qi,DPi)
        V4,D3,Bv,R4,DR,Xtet,Ytet = volute(V1,Ws,rot,DPi,Qi,B1,D2,alpha2)
        #Function of Blade design
        Rb = blades(D1,D2,beta1)
        if QL == 0:
            #First input data not right
            print("\n    * The input data: Q = {:1.4f} mc/s - DP = {:1.4f} Pa".format(Qi,DPi))
            print("    => Correction of input data with hydraulic, leakage and power losses")
            #Correction with hydraulic, leakage and power losses
            #Leakage losses
            Cd = 0.65 #The disharge coefficient
            cl = 0.01 * D1 #clearance (1% D1)
            QL = Cd * (pi * D1) * cl * num.sqrt((4*abs(DPi)/3)/rot)
            print("       - Flow rate correction: {:3.5f} cm/s".format(QL))
            #Suction pressure loss - ki - loss factor 0.1
            ki = 0.1
            dpsuct = 0.5 * rot * ki * Vduct**2
            print("       - Suction pressure loss: {:3.5f} Pa".format(dpsuct))
            #Impeller pressure loss - kii - loss factor 0.2 - 0.3
            kii = 0.25
            dpimp = 0.5 * rot * kii * (W1-W2)**2
            print("       - Impeller pressure loss: {:3.5f} Pa".format(dpimp))
            #Volute pressure loss - kiii - loss factor 0.2
            kiii = 0.3
            dpvol = 0.5 * rot * kiii * (V2-V4)**2
            print("       - Volute pressure loss: {:3.5f} Pa".format(dpvol))

            ##Discharge and pressure corrections + Efficiency estimation
            Qii = Qi + QL
            DPii = abs(DPi) + dpsuct + dpimp + dpvol
            etaHyd = abs(DPi)/DPii
            etaVol = Qi/Qii
            etaTot = etaHyd * etaVol

            #corrected values of the inputa data
            phi2 = (W2/U2)
            slip = 1 - (1.98/Nb)/(1-phi2/num.tan(beta2)) #Stanitz Formula complete formula
            Qi = Qii
            DPi = DPii
        else:
            print("    * Exact input data: Q = {:1.4f} mc/s - DP = {:1.4f} Pa".format(Qi,DPi))
            print("    * Impeller diameters: D1 = {:1.4f} m - D2 = {:1.4f} m".format(D1,D2))
            ef ="    * Estimated Efficiency: etaHydra = {:8.6f} - etaVolum = {:8.6f} - etaTot = {:8.6f}"
            print(ef.format(etaHyd,etaVol,etaTot))
            Qf = Qi
            DPf = DPi
            corr = True

    #Estimation of the shaft diameter Dshaft
    #Disk friction loss - f = 0.005 friction factor for mild steel sheet
    f = 0.005
    Tdf = pi * rot * f * ((2*U2/D2)**2) * ((D2/2)**5)/5
    Pdf = w*Tdf
    Pideal = (DPf*Qf)/etaTot + Pdf
    Tideal = Pideal/w
    safe = 4 #safety coefficient
    tau = 343e+5
    Dshaft = (16*Tideal*safe/(pi * tau))**(1/3)

    #output of the Design function
    Ainlet =  [D1,Dduct,Lduct,U1,V1,Vduct,W1,beta1,B1,A1,Aduct]
    Aoutlet = [U2,D2,V2,V2m,W2,alpha2,beta2,B2,A2,Pow,Ws,slipf]
    Avolute = [V4,D3,Bv,R4,DR,Xtet,Ytet]
    Ac      = [Qf,DPf,etaHyd,etaVol,etaTot,Pideal,Tideal,Dshaft]

    return [Ainlet, Aoutlet, Avolute,Ac,Rb]

#############################
def thisIsMyDesign(nameF,Qf,DPf,etaHyd,etaVol,etaTot,Pideal,Tideal,Dshaft,D1,Dduct,Lduct,B1,Nb,D2,beta2,B2,D3,Bv,R4,DR,Xtet,Ytet):
    #Writing on a file
    name = "Design_" + nameF + ".txt"
    data = open(name,'w')
    data.write("###########################################################\n")
    data.write("###########################################################\n")
    data.write("#########                                         #########\n")
    data.write("#########       Design Parameters                 #########\n")
    data.write("#########          of a Centrifugal Fan           #########\n")
    data.write("#########                                         #########\n")
    data.write("#########     NAME: {:16s} ver.1.2      #########\n".format(nameF))
    data.write("###########################################################\n")
    data.write("###########################################################\n")
    data.write('Qf  = {:3.5f} cm/s\n'.format(Qf))
    data.write("DPs = {:3.5f} Pa\n".format(DPf))
    data.write("Pid = {:3.5f} W\n".format(Pideal))
    data.write("D1  = {:3.5f} m - Inner diameter\n".format(D1))
    data.write("D2  = {:3.5f} m - Outer diameter\n".format(D2))
    data.write("B1  = {:3.5f} m - inner width\n".format(B1))
    data.write("B2  = {:3.5f} m - outer width\n".format(B2))
    data.write("Nb  = {:3.5f} - Number of blades\n".format(Nb))
    data.write("etaT  = {:3.5f} - Total Efficiency\n".format(etaTot))
    data.write("etaH  = {:3.5f} - Hydraulic Efficiency\n".format(etaHyd))
    data.write("etaV  = {:3.5f} - Volumetric Efficiency\n".format(etaVol))
    data.close()



def inPutParam():
    rad = (2*pi/60) #coversion factor in rad/s (omega)
    mil = 1000 #conversion factor in mm
    nameF = input("Name of the fan: ") #Assign a name to the Fan
    t = float(input("- Select the air temperature t(C): ")) #Set an operating temperature
    qo = float(input("- Select the air flow rate Q(m^3/s): ")) #Set a flow rate
    pin = float(input("- Select the inlet static pressure Pin(Pa): "))
    pout = float(input("- Select the outlet static pressure Pout(Pa): "))
    dPo = abs(pin - pout) # Static Pressure Drop
    nb = float(input("- Select the number of blades Nb(-): "))
    th = float(input("- Select the blade thickness th(mm): "))
    w = float(input("- Select the rotational speed w(rpm): "))

    th = th/mil # thickness in meter
    w = w*rad   # rot. speed in rad/sec
    T = t + 273 # Temperature in Kelvin
    return nameF,qo,pin,dPo,nb,th,w,T

##############################################################################
###   Main    #####
#Design input parameter
nameF,Qo,Pin,DPo,Nb,th,w,T = inPutParam()
Qi = Qo
DPi = DPo
#Design Function
A = design(nameF, Qi, Pin, DPi, Nb, th, w, T)
#Design output
D1,Dduct,Lduct,U1,V1,Vduct,W1,beta1,B1,A1,Aduct = A[0]
U2,D2,V2,V2m,W2,alpha2,beta2,B2,A2,Pow,Ws,slipf = A[1]
V4,D3,Bv,R4,DR,Xtet,Ytet = A[2]
Qf,DPf,etaHyd,etaVol,etaTot,Pideal,Tideal,Dshaft = A[3]
Rb = A[4]

#WriteFunction
thisIsMyDesign(nameF,Qf,DPf,etaHyd,etaVol,etaTot,Pideal,Tideal,Dshaft,D1,Dduct,Lduct,B1,Nb,D2,beta2,B2,D3,Bv,R4,DR,Xtet,Ytet)
input("Press enter to exit!")

#A[0] = [D1,Dduct,U1,V1,V1m,W1,beta1,B1,A1,Aduct]
#A[1] = [U2,D2,V2,V2m,W2,alpha2,beta2,B2,A2,Pow,Ws]
##A[2] = [V4,D3,Bv,R4,DR,Xtet,Ytet]
##A[3] = [Qf,DPf,etaHyd,etaVol,etaTot,Pideal,Tideal,Dshaft]
##A[4] = Rb
