# -*- coding: utf-8 -*-
"""
Created on Wed March 14 09:21:55 2018

Last Update: 26/09/2018 
ver. 1
This version allow to estimate the Fan geometry taking into account the 
compressibility of the flow. The blade has the same width.

ver. 2
The affinity laws are added
 
To-Do list:
- Compressible effect on the outer diameter
            
Note1 >> For Boundary layer use the tool BoundaryLayer.py !!

@author: Vincenzo Sammartano
@buchermunicipal
"""

### Libraries
import numpy as num
from math import pi
#import matplotlib.pyplot as plt
#import sys
#import shutil

### Functions
def impellerIn(Qi,rot,Tatm,uvRatio):
    #inlet Duct and impeller inlet design
    #Duct and impeller parameters
    teta = num.radians(30) # duct contraction
    ind = 0.2  #induction improvement
    Q1 = 0
    iq = 0
    while (Q1 == 0):
        Q1 = Qi
        #Duct geometry and Impeller eye geometry
        D1 = num.power((uvRatio * 8/pi) * (Q1/w),1/3)
        Dduct = ( 1+ind ) * D1
        
        A1 = 0.25 * pi * D1**2
        Aduct = 0.25 * pi * Dduct**2
        Lduct = 0.1 * D1 / (2 * num.tan(teta))
    
        #Velocity triangle
        U1 = w * D1/2
        V1 = U1/uvRatio
        
        Vduct = 4*Q1/(pi*(Dduct**2))
        W1 = num.power(U1**2+V1**2,0.5)
        beta1 = num.rad2deg(num.arctan(1/uvRatio))
        #impeller width at inlet section
    
        corr1 = 1 # correction factor for blade thickness
        B1 = (Q1/V1)/(pi*D1 - Nb*th)/corr1
        
        Pin = 101325 - 9.806 * rot * (V1**2)/(2*9.806) #Pressure inlet
        epIn = Pin / 101325 # pressure ratio
        Tin = (Tatm * epIn**0.286)
        Rf = 287.058
        rot1 = Pin/(Rf*Tin)
        Qi = Q1 * rot/rot1
        iq += 1

        if iq == 1: Q1 = 0
        
    return [D1,Dduct,Lduct,U1,V1,Vduct,W1,beta1,B1,A1,Aduct,Qi,rot1]

def impellerOut(DPi,Qi,rot,slip,B1):
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

def volute(c2,V1,Ws,rot,DPi,Qi,B1,D2,alpha2):
    clv = c2 * D2 # clearance = 1.5% D2
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

def design(nameF,Qi,Pin,DPi,Nb,th,w,T,c1,c2,uvRatio):
    Rf = 287.058 # Universal Constant of Gases [J/(Kg K)]
    pAtm = 101325 # [Pa] atmospheric pressure
    slip = 1-(1.98/Nb) #Stanitz
    #Taking compressibility effect into consideration
    rot = pAtm/(Rf*T) #air density at T celsius
    corr = False
    QL = 0 #Flow rate

    while corr == False:
        #design with the input data Q & DP
        D1,Dduct,Lduct,U1,V1,Vduct,W1,beta1,B1,A1,Aduct,Q1,rot1 = impellerIn(Qi,rot,T,uvRatio)
        Qi = Q1 #corrected value of the flow area due to the compressibility effect
        rot = rot1 #the new density value due to the compressibility effect
        U2,D2,V2,V2m,W2,alpha2,beta2,B2,A2,Pow,Ws,slipf = impellerOut(DPi,Qi,rot,slip,B1)
        #function to estimate the volute casing geometry - F(Qi,DPi)
        V4,D3,Bv,R4,DR,Xtet,Ytet = volute(c2,V1,Ws,rot,DPi,Qi,B1,D2,alpha2)
        #Function of Blade design
        Rb = blades(D1,D2,beta1)
        
        if QL == 0:
            #First input data not right
            print("\n ... correcting the data input with hydraulic, leakage and power losses:")
            #Correction with hydraulic, leakage and power losses
            #Leakage losses
            Cd = 0.65 #The discharge coefficient
            cl = c1 * D1 #clearance (1% D1)
            QL = Cd * (pi * D1) * cl * num.sqrt((4*abs(DPi)/3)/rot)
            print("\n       - Flow rate correction: {:3.5f} cm/s".format(QL))
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
            if slip < 1:
                slip = 1 - (1.98/Nb)
            
            Qi = Qii
            DPi = DPii
        else:
            print("\n- Impeller diameters: D1 = {:1.4f} m - D2 = {:1.4f} m".format(D1,D2))
            print("\n- Impeller width: B1 = {:1.4f} m - B2 = {:1.4f} m".format(B1,B2))
            print("\n- Rotational Speed: omega = {:4.0f} rpm".format(w/rad))
            ef = "\n- Estimated Efficiency:\n\tetaHydra = {:8.6f}\n\tetaVolum = {:8.6f}\n\tetaTot = {:8.6f}\n"
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
    Ainlet =  [D1,Dduct,Lduct,U1,V1,Vduct,W1,beta1,B1,A1,Aduct,rot1] 
    Aoutlet = [U2,D2,V2,V2m,W2,alpha2,beta2,B2,A2,Pow,Ws,slipf]
    Avolute = [V4,D3,Bv,R4,DR,Xtet,Ytet]
    Ac      = [Qf,DPf,etaHyd,etaVol,etaTot,Pideal,Tideal,Dshaft]

    return [Ainlet, Aoutlet, Avolute,Ac,Rb]

def affinityLwas(Q,Dp,w,D):
    text = """
    ~~~~~~~ Affinity laws of the Fan ~~~~~~~~
      1)    Q1/Q2 = (w1/w2) * (d1/d2)^3
    
      2)    dp1/dp2 = (w1/w2)^2 *(d1/d2)^2
    
      3)    P1/P2 = (w1/w2)^3 * (d1/d2)^5
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    """
    print(text)

    q1 = Q
    d1 = D
    w1 = w/rad #in rpm
    dp1 = Dp

    affin = 999
    while (affin not in range(1,3)):
        affin = int(input("\t1. Changing impeller Velocity\n\t2. Changing the Impeller Diameter\n\t- Please select [1-2]: "))
        if affin == 1:
            print("\n***** Changing impeller Velocity *****")
            w2 = float(input("- Select the new impeller speed omega(rpm) = "))
            d2 = d1
            q2 = q1 * ((w1/w2)**-1) * ((d1/d2)**-3)
            dp2 = dp1 * ((w1/w2)**-2) * ((d1/d2)**-2)
            P2 = q2*dp2/mil #in KWatt
        elif affin == 2:
            print("\n***** Changing impeller Diameter *****")
            d2 = float(input("- Select the new impeller diameter(mm) = "))
            d2 = d2/mil
            w2 = w1
            q2 = q1 * ((d1/d2)**-3)
            dp2 = dp1 * ((d1/d2)**-2)
            P2 = q2*dp2/mil #in KWatt
        else:
            print("\n\t ... your selection is wrong!!!\n\t ... please select [1-2]")

    print('\n- New impeller speed omega\t= {:8.1f} rpm'.format(w2))
    print('- New impeller diameter\t= {:6.1f} mm'.format(d2*mil))
    print('- New flow rate Q\t= {:8.4f} m^3/s'.format(q2))
    print('- New pressure rise DP\t= {:8.2f} Pa'.format(dp2))
    print('- New generated Power Pow\t= {:8.2f} kW\n'.format(P2))

    return q2,dp2,P2,w2,d2 #Qi_new,DPi_new,P_new,w_new,D1_new


#############################
def thisIsMyDesign(nameF,Qf,DPf,etaHyd,etaVol,etaTot,Pideal,Tideal,Dshaft,
                   beta1,D1,Dduct,Lduct,B1,Nb,D2,beta2,B2,D3,Bv,R4,DR,Xtet,Ytet,
                   Qi_new,DPi_new,P_new,w_new,D1_new,aff):

    #Writing on a file
    name = "Design_" + nameF + ".txt"
    data = open(name,'w')
    data.write("###########################################################\n")
    data.write("###########################################################\n")
    data.write("#########                                         #########\n")
    data.write("#########       Design Parameters                 #########\n")
    data.write("#########          of a Centrifugal Fan           #########\n")
    data.write("#########                                         #########\n")
    data.write("#########     NAME: {:16s} ver.{}      #########\n".format(nameF,ver))
    data.write("###########################################################\n")
    data.write("###########################################################\n\n")
    data.write('  * T  = {:> 8.1f} C - Air Flow temperature\n'.format(t))
    data.write('  * Q  = {:> 8.4f} m^3/s - Air Flow rate\n'.format(Qi))
    data.write("  * DP = {:> 8.2f} Pa - Pressure rise\n".format(DPi))
    data.write("  * Pgen = {:> 8.2f} KW - Generated Power\n".format(Qi*DPi/mil))
    data.write("  * Speed = {:> 8.0f} rpm - Rotational speed\n\n".format(w/rad))

    data.write("  * Dduct  = {:> 8.2f} mm - Duct diameter\n".format(Dduct*mil))
    data.write("  * Lduct  = {:> 8.2f} mm - Duct lenght\n\n".format(Lduct*mil))

    data.write("  * D1\t= {:> 8.2f} mm - Inner diameter\n".format(D1*mil))
    data.write("  * D2\t= {:> 8.2f} mm - Outer diameter\n".format(D2*mil))
    data.write("  * B1\t= {:> 8.2f} mm - inner width\n".format(B1*mil))
    data.write("  * B2\t= {:> 8.2f} mm - outer width\n".format(B2*mil))
    data.write("  * beta1 = {:> 8.1f} deg - attack angle (U1^W1)\n".format(beta1))
    data.write("  * beta2 = {:> 8.1f} deg - leaving angle (U2^W2)\n\n".format(beta2))

    data.write("  * Cl1\t= {:> 8.1f}% {:3.1f} = {:3.1f} mm - Inlet cleareances\n".format(c1*cent,D1*mil,c1*D1*mil))
    data.write("  * Cl2\t= {:> 8.1f}% {:3.1f} = {:3.1f} mm - Outlet cleareances\n".format(c2*cent,D2*mil,c2*D2*mil))
    data.write("  * Nb\t= {:> 8.0f} - Number of blades\n".format(Nb))
    data.write("  * th\t= {:> 8.1f} mm - Blades thickness \n\n".format(th*mil))
    
    data.write("  * etaH  = {:> 8.2f}% - Flow Efficiency\n".format(etaHyd*cent))
    data.write("  * etaV  = {:> 8.2f}% - Volumetric Efficiency\n\n".format(etaVol*cent))
    data.write("  * etaT  = {:> 8.2f}% - Static Efficiency\n".format(etaTot*cent))

    if aff == True:
        data.write('\n\n---------------------  Affinity law ---------------------\n')
        #Qi_new,DPi_new,P_new,w_new,D1_new
        data.write('\n  * New impeller speed omega\t= {:> 8.0f} rpm\n'.format(w_new))
        data.write('\n  * New impeller diameter\t= {:> 8.2f} mm\n'.format(D2_new*mil))
        data.write('\n  * New flow rate Q\t= {:> 8.4f} m^3/s\n'.format(Qi_new))
        data.write('\n  * New pressure rise DP\t= {:> 8.1f} Pa\n'.format(DPi_new))
        data.write('\n  * New generated Power Pow\t= {:> 8.2f} kW\n'.format(P_new))
        data.write('\n-------------------------------------------------------------')
    data.close()

def inPutParam():
    print("###################################################################")
    print("#########         Centrifugal Fan Design ver.{}           ########".format(ver))
    print("###################################################################\n")

    nameF = input("Name of the fan: ") #Assign a name to the Fan
    t = float(input("\n- Select the air temperature t(C): ")) #Set an operating temperature
    qo = float(input("\n- Select the air flow rate Q(m^3/s): ")) #Set a flow rate
    pin = float(input("\n- Select the inlet static pressure Pin(Pa): "))
    pout = float(input("\n- Select the outlet static pressure Pout(Pa): "))
    dPo = abs(pin - pout) # Static Pressure Drop
    nb = float(input("\n- Select the number of blades Nb(-): "))
    th = float(input("\n- Select the blade thickness th(mm): "))
    w = float(input("\n- Select the rotational speed w(rpm): "))

    c1 = float(input("\n- Select the cleareances % of the Din  c1 (1% D1)   =  "))
    c2 = float(input("\n- Select the cleareances % of the Dout c2 (1.5% D2) =  "))
    alfAtt = float(input("\n- Select the attack angle - beta1 (35-45 deg):  "))
    uvRatio = 1/num.tan(num.deg2rad(alfAtt))

    th = th/mil # thickness in meter
    w = w*rad   # rot. speed in rad/sec
    T = t + 273 # Temperature in Kelvin
    return nameF,qo,pin,dPo,nb,th,w,t,T,c1,c2,uvRatio

##############################################################################
###   Main    #####
#Design input parameter
ver = "2.0" #version of the code
cent = 100
mil = 1000 #conversion factor in mm
rad = (2*pi/60) #coversion factor in rad/s (omega)
nameF,Qo,Pin,DPo,Nb,th,w,t,T,c1,c2,uvRatio = inPutParam()
Qi = Qo
DPi = DPo

#Design Function
A = design(nameF, Qi, Pin, DPi, Nb, th, w, T,c1,c2,uvRatio)

#Design output
D1,Dduct,Lduct,U1,V1,Vduct,W1,beta1,B1,A1,Aduct,rot1 = A[0]
U2,D2,V2,V2m,W2,alpha2,beta2,B2,A2,Pow,Ws,slipf = A[1]
V4,D3,Bv,R4,DR,Xtet,Ytet = A[2]
Qf,DPf,etaHyd,etaVol,etaTot,Pideal,Tideal,Dshaft = A[3]
Rb = A[4]

#Affinity laws
aff = input("Would like to use the affinity laws? [y/n] ... ")
if aff in ['Y','y','Yes','yes','YES']:
    aff = True
    #Affinity laws
    Qi_new,DPi_new,P_new,w_new,D2_new = affinityLwas(Qi,DPi,w,D2)
else:
    aff = False
    [Qi_new,DPi_new,P_new,w_new,D2_new] = [0,0,0,0,0]
    
#WriteFunction
thisIsMyDesign(nameF,Qf,DPf,etaHyd,etaVol,etaTot,Pideal,Tideal,Dshaft,beta1,
               D1,Dduct,Lduct,B1,Nb,D2,beta2,B2,D3,Bv,R4,DR,Xtet,Ytet,
               Qi_new,DPi_new,P_new,w_new,D2_new,aff)

input("Press enter to exit!")

