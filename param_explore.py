##COMPENSATION POINT FUNCTION

import sys
import numpy as np                         # Numpy for some math stuff
from scipy.integrate import odeint         # the odeint function is what we'll use to solve our system of ODEs
import pandas as pd
import math
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
import time


def final_fluxes(y0, T1, T2, T3, Pmc, Pcp, Pup, Puc, Po, CO2m, O2m, HCm, CA_conc, Vc, Kc, Ko, k2, k_2, ShCO2_forward, ShCO2_reverse, khcat, KCO2, KHC, Vmax_BicA, Km_BicA, kin_Sco, Keq_cyt, Keq_cpl, Pmb, T4, pH_chloroplast, pH_cytosol, q10_PlipCO2, q10_PlipHCO3, PlipHCO3, pumpCost, PlipCO2, radCell, protonc, protonp, SA_Cplast, Nmem, Pmo):
    CO2c,CO2b,HCc,CO2p,HCp,O2p,O2b = y0

    # CO2 diffusing in from outside
    
    V18 = Pmb*(CO2m - CO2b)*1e-6
    
    V19 = Pmo*(O2m - O2b)*1e-6
    
    V1 = Pmc*(CO2b - CO2c)*1e-6
    
    # Spontaneous hydration (forward reverse) of CO2 in cytosol
    
    V2 = CO2c * k2 * 1e-6 * T3 
    V3 = HCc * protonc * k_2 * 1e-6 * T3 
    
    # CA mediated interconversion of CO2 and HCO3 (cytosol)
    
    V4 = (khcat*(CA_conc)*(((CO2c)-(((HCc)*(protonc)/(1.071519305237607))))))/((KCO2)+((KCO2/KHC)*(HCc))+(CO2c)) * T3   
    
    # Spontaneous movement of CO2 through chloroplast membrane
    
    V5 = Pcp * (CO2c - CO2p) * 1e-6 
    
    # Spontaneous movement of HCO3 through plasmalemma
    
    V6 = Puc * (HCm - HCc) * 1e-6  
    
    # Spontaneous movement of HCO3- through chloroplast membrane

    V7 = Pup * (HCc - HCp) * 1e-6 
    
    # Active transport of HCO3- through chloroplast membrane
    
    V8 = ((Vmax_BicA * HCc)/(HCc + Km_BicA)) * SA_Cplast
    
    # Spontaneous hydration (forward reverse) of CO2 in stroma
    
    V9 = CO2p * k2 * 1e-6 * T1 
    V10 = HCp * protonp * k_2 * 1e-6 * T1 
    
    # CA mediated interconversion of CO2 and HCO3 (stroma)
    
    V11 = (khcat*(CA_conc)*(((CO2p)-(((HCp)*(protonp)/(1.071519305237607))))))/((KCO2)+((KCO2/KHC)*(HCp))+(CO2p)) * T1    #CA-catalyzed interconversion of bicarb and CO2 in stroma 
    
    # CO2 fixation by RuBisCO
    
    V12 = (Vc*CO2p)/(CO2p+(Kc*((1+(O2p/Ko))))) 
    
    # Drawdown of O2 by RuBisCO
    
    vc_vo = kin_Sco * (CO2p/O2p) 
    V13 = V12 / vc_vo 
  
    # Evolution of CO2 in cytoplasm as a result of photorespiration
    
    V14 = 0.5 * V13 
    
    # Diffusion of O2 out of chloroplast
    
    V15 = Po*(O2p - O2b) * 1e-6 
    
    # Evolution of O2
    
    V16 = V12  
    
    # Despiration in the light (RL)
    
    V17 = 4.95e-19 * (T3/3.730641276137879e-15) #mol CO2 s^-1 cell^-1 (magnitude of respiration in the light) assuming 1um cell radius
    
    # Calculate ATP cost
    
    ##from photosynthesis:
    ATP_PS = (V12*1e18)*3
    ##from photorespiration:
    ATP_PR = (V13*1e18)*3.5
    ##from bicarbonate pumping into stroma:
    ATP_pump = pumpCost * (V8*1e18)
    ##total cost of operating CCM
    ATP_total = ATP_PS + ATP_PR + (ATP_pump*Nmem)
    
    #Calculate NADPH cost
    ##from photosynthesis and photorespiration:
    NADPH_total = 2*(V12+V13)*1e18
    
    #calculate net assimilation 
    assim_net = V12-V14-V17
    
    # Calculating changes in concentration based off these rates
    
    dCO2cdt = 1e6*((V1 + V3 + V14 + V17 - V2 - V4 - V5)/T3)
    dCO2bdt = 1e6*((V18 - V1)/T4)
    dHCO3cdt = 1e6*((V2 + V4 + V6 - V3 - V7 - V8)/T3)
    dCO2pdt = 1e6*((V5 + V10 - V11 - V12 - V9)/T1)
    dHCO3pdt = 1e6*((V7 + V8 + V9 + V11 - V10)/T1)
    dO2pdt = 1e6*((V16 - V13 - V15)/T1)
    dO2bdt = 1e6*((V19 + V15)/T4)
    
    return assim_net, 1/vc_vo, ATP_total, NADPH_total, V5/V12, ATP_total/(assim_net*1e18)

t = np.linspace(0,1,100000)


CO2b = 0.1
CO2p = 0.1            # uM 
CO2c = 0.0            # uM
O2p = 0.1             # uM
HCc = 0.001           # uM
HCp = 0.001           # uM 
O2b = 0.1             # uM

y0 = [CO2c,CO2b,HCc,CO2p,HCp,O2p,O2b]



def dudx(y0, t, T1, T2, T3, Pmc, Pcp, Pup, Puc, Po, CO2m, O2m, HCm, CA_conc, Vc, Kc, Ko, k2, k_2, ShCO2_forward, ShCO2_reverse, khcat, KCO2, KHC, Vmax_BicA, Km_BicA, kin_Sco, Keq_cyt, Keq_cpl, Pmb, T4, pH_chloroplast, pH_cytosol, q10_PlipCO2, q10_PlipHCO3, PlipHCO3, pumpCost, PlipCO2, radCell, protonc, protonp, SA_Cplast, Nmem, Pmo):
    CO2c,CO2b,HCc,CO2p,HCp,O2p,O2b = y0

    # CO2 diffusing in from outside
    
    V18 = Pmb*(CO2m - CO2b)*1e-6
    
    V19 = Pmo*(O2m - O2b)*1e-6
    
    V1 = Pmc*(CO2b - CO2c)*1e-6
    
    # Spontaneous hydration (forward reverse) of CO2 in cytosol
    
    V2 = CO2c * k2 * 1e-6 * T3 
    V3 = HCc * protonc * k_2 * 1e-6 * T3 
    
    # CA mediated interconversion of CO2 and HCO3 (cytosol)
    
    V4 = (khcat*(CA_conc)*(((CO2c)-(((HCc)*(protonc)/(1.071519305237607))))))/((KCO2)+((KCO2/KHC)*(HCc))+(CO2c)) * T3   
    
    # Spontaneous movement of CO2 through chloroplast membrane
    
    V5 = Pcp * (CO2c - CO2p) * 1e-6 
    
    # Spontaneous movement of HCO3 through plasmalemma
    
    V6 = Puc * (HCm - HCc) * 1e-6  
    
    # Spontaneous movement of HCO3- through chloroplast membrane

    V7 = Pup * (HCc - HCp) * 1e-6 
    
    # Active transport of HCO3- through chloroplast membrane
    
    V8 = ((Vmax_BicA * HCc)/(HCc + Km_BicA)) * SA_Cplast
    
    # Spontaneous hydration (forward reverse) of CO2 in stroma
    
    V9 = CO2p * k2 * 1e-6 * T1 
    V10 = HCp * protonp * k_2 * 1e-6 * T1 
    
    # CA mediated interconversion of CO2 and HCO3 (stroma)
    
    V11 = (khcat*(CA_conc)*(((CO2p)-(((HCp)*(protonp)/(1.071519305237607))))))/((KCO2)+((KCO2/KHC)*(HCp))+(CO2p)) * T1    #CA-catalyzed interconversion of bicarb and CO2 in stroma 
    
    # CO2 fixation by RuBisCO
    
    V12 = (Vc*CO2p)/(CO2p+(Kc*((1+(O2p/Ko))))) 
    
    # Drawdown of O2 by RuBisCO
    
    vc_vo = kin_Sco * (CO2p/O2p) 
    V13 = V12 / vc_vo 
  
    # Evolution of CO2 in cytoplasm as a result of photorespiration
    
    V14 = 0.5 * V13 
    
    # Diffusion of O2 out of chloroplast
    
    V15 = Po*(O2p - O2b) * 1e-6 
    
    # Evolution of O2
    
    V16 = V12  
    
    # Despiration in the light (RL)
    
    V17 = 4.95e-19 * (T3/3.730641276137879e-15) #mol CO2 s^-1 cell^-1 (magnitude of respiration in the light) assuming 1um cell radius
    
    # Calculate ATP cost
    
    ##from photosynthesis:
    ATP_PS = (V12*1e18)*3
    ##from photorespiration:
    ATP_PR = (V13*1e18)*3.5
    ##from bicarbonate pumping into stroma:
    ATP_pump = pumpCost * (V8*1e18)
    ##total cost of operating CCM
    ATP_total = ATP_PS + ATP_PR + (ATP_pump*Nmem)
    
    #Calculate NADPH cost
    ##from photosynthesis and photorespiration:
    NADPH_total = 2*(V12+V13)*1e18
    
    #calculate net assimilation 
    assim_net = V12-V14-V17
    
    # Calculating changes in concentration based off these rates
    
    dCO2cdt = 1e6*((V1 + V3 + V14 + V17 - V2 - V4 - V5)/T3)
    dCO2bdt = 1e6*((V18 - V1)/T4)
    dHCO3cdt = 1e6*((V2 + V4 + V6 - V3 - V7 - V8)/T3)
    dCO2pdt = 1e6*((V5 + V10 - V11 - V12 - V9)/T1)
    dHCO3pdt = 1e6*((V7 + V8 + V9 + V11 - V10)/T1)
    dO2pdt = 1e6*((V16 - V13 - V15)/T1)
    dO2bdt = 1e6*((V19 + V15)/T4)
    
    return dCO2cdt, dCO2bdt, dHCO3cdt, dCO2pdt, dHCO3pdt, dO2pdt, dO2bdt

def compensation_calculation(args):

    steady_state_check = []
    CO2 = np.logspace(-4,3,num=100,endpoint=True)
    net_assimilation_values = []

    for x in range(len(CO2)):
            # Check here 
            solution = odeint(dudx, y0, t, mxstep=5000, args=(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], CO2[x], args[9], args[10], args[11], args[12], args[13], args[14], args[15], args[16], args[17], args[18], args[19], args[20], args[21], args[22], args[23], args[24], args[25], args[26], args[27], args[28], args[29], args[30], args[31], args[32], args[33], args[34], args[35], args[36], args[37], args[38], args[39], args[40], args[41]))
            steady_state_concs = [solution[-1][0],solution[-1][1],solution[-1][2],solution[-1][3],solution[-1][4],solution[-1][5],solution[-1][6]]
            net_assim_final, vo_vc_final, ATP_final, NADPH_final, leak_final, ATP_per_CO2_final = final_fluxes(steady_state_concs, args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10], args[11], args[12], args[13], args[14], args[15], args[16], args[17], args[18], args[19], args[20], args[21], args[22], args[23], args[24], args[25], args[26], args[27], args[28], args[29], args[30], args[31], args[32], args[33], args[34], args[35], args[36], args[37], args[38], args[39], args[40], args[41])
            net_assimilation_values.append(net_assim_final)
            if solution[-1][0] <= solution[-2][0]*1.0001 and solution[-1][0] >= solution[-2][0]*0.9999 and solution[-1][1] <= solution[-2][1]*1.0001 and solution[-1][1] >= solution[-2][1]*0.9999 and solution[-1][2] <= solution[-2][2]*1.0001 and solution[-1][2] >= solution[-2][2]*0.9999 and solution[-1][3] <= solution[-2][3]*1.0001 and solution[-1][3] >= solution[-2][3]*0.9999 and solution[-1][4] <= solution[-2][4]*1.0001 and solution[-1][4] >= solution[-2][4]*0.9999 and solution[-1][5] <= solution[-2][5]*1.0001 and solution[-1][5] >= solution[-2][5]*0.9999 and solution[-1][6] <= solution[-2][6]*1.0001 and solution[-1][6] >= solution[-2][6]*0.9999:
                steady_state_check.append(0)
                continue
            else:
                print('Warning! Steady state not reached')
                steady_state_check.append(1)
    spline = CubicSpline(CO2,net_assimilation_values)
    root = spline.roots()[-1]
    
    # Now calculate at CO2m
    
    solution = odeint(dudx, y0, t, mxstep=5000, args=(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10], args[11], args[12], args[13], args[14], args[15], args[16], args[17], args[18], args[19], args[20], args[21], args[22], args[23], args[24], args[25], args[26], args[27], args[28], args[29], args[30], args[31], args[32], args[33], args[34], args[35], args[36], args[37], args[38], args[39], args[40], args[41]))
    if solution[-1][0] <= solution[-2][0]*1.0001 and solution[-1][0] >= solution[-2][0]*0.9999 and solution[-1][1] <= solution[-2][1]*1.0001 and solution[-1][1] >= solution[-2][1]*0.9999 and solution[-1][2] <= solution[-2][2]*1.0001 and solution[-1][2] >= solution[-2][2]*0.9999 and solution[-1][3] <= solution[-2][3]*1.0001 and solution[-1][3] >= solution[-2][3]*0.9999 and solution[-1][4] <= solution[-2][4]*1.0001 and solution[-1][4] >= solution[-2][4]*0.9999 and solution[-1][5] <= solution[-2][5]*1.0001 and solution[-1][5] >= solution[-2][5]*0.9999 and solution[-1][6] <= solution[-2][6]*1.0001 and solution[-1][6] >= solution[-2][6]*0.9999:
        steady_state_check.append(0)
    else:
        print('Warning! Steady state not reached')
        steady_state_check.append(1)
    
    if sum(steady_state_check) == 0: 
        steady_state_result = 'Pass'
    else: 
        steady_state_result = 'Fail'
        
    # Calculating relevant fluxes
    
    steady_state_concs = [solution[-1][0],solution[-1][1],solution[-1][2],solution[-1][3],solution[-1][4],solution[-1][5],solution[-1][6]]
    net_assim_final, vo_vc_final, ATP_final, NADPH_final, leak_final, ATP_per_CO2 = final_fluxes(steady_state_concs, args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8], args[9], args[10], args[11], args[12], args[13], args[14], args[15], args[16], args[17], args[18], args[19], args[20], args[21], args[22], args[23], args[24], args[25], args[26], args[27], args[28], args[29], args[30], args[31], args[32], args[33], args[34], args[35], args[36], args[37], args[38], args[39], args[40], args[41])
    
    return root, solution[-1][0], solution[-1][1], solution[-1][2], solution[-1][3], solution[-1][4], solution [-1][5], solution[-1][6], net_assim_final*1e18, vo_vc_final, ATP_final, NADPH_final, leak_final, ATP_per_CO2, steady_state_result

def parameter_calculation(CA_conc = 0.27, khcat = 0.3e6, KHC = 34, KCO2 = 1.5, pumpCost = 1, Vmax_BicA = 185e-6, PlipHCO3 = 2e-9, radCell = 2.0, q10_PlipHCO3 = 1, q10_PlipCO2 = 1, pH_cytosol = 7, pH_chloroplast = 8, Km_BicA = (0.217/1000)*1e6, generic_Q10=2,Nmem=5,temperature=45,PlipCO2=3.5e3,BicA_Vmax_factor=1,radBoundary=10,kin_Sco = 238,Vc = 0.825e-16 * 4.5,Kc = 24.9, Ko = 479):
    
    # -- calculate cell size parameters  -- 
        
    radCell = radCell #um (the radius of the cell)
    radBoundary = radCell*2 
    radCplast = radCell/2 #um (the radius of the chloroplast)
    radPyr = radCell/4  #um (the radius of the clear space in middle of chloroplast)
    SA_Pyr = 4*np.pi*(radPyr**2) #um^2 (the surface area of the clear space in middle of chloroplast)
    SA_Cplast = 4*np.pi*(radCplast**2) #um^2 (the surface area of the chloroplast)
    SA_Cell = 4*np.pi*(radCell**2) #um^2 (the surface area of the cell)
    SA_Boundary = 4*np.pi*(radBoundary**2)
    T1 = (4/3)*np.pi*(radPyr**3) / 1e15 # um^3, converted to L (stroma space volume)
    T2 = (4/3)*np.pi*(radCplast**3) / 1e15 - T1 # um^3, converted to L (chloroplast volume without the inner stroma space)
    T3 = (4/3)*np.pi*(radCell**3)  / 1e15 - T2  # um^3, converted to L (cytoplasm volume)
    T4 = (4/3)*np.pi*(radBoundary**3) / 1e15 - T3
    LengthBoundary = radBoundary - radCell
    
    # -- input temperature adjustment parameters --

    temperature = temperature #C
    q10_PlipCO2 = q10_PlipCO2
    q10_PlipHCO3 = q10_PlipHCO3
    q10_Ka_HCtransport = 2
    q10_Vc = 2.21 #vonCaemmerer book (Vcmax)
    q10_Kc = 2.24 #vonCaemmerer book
    q10_Ko = 1.63 #vonCaemmerer book
    q10_Kba = 2
    q10_k2 = 2 
    q10_k_2 = 2  
    q10_ShCO2_forward = 2
    q10_ShCO2_reverse = 2 
    q10_khcat = 2 
    q10_KCO2 = 2 
    q10_KHC = 2 
    q10_Vmax_BicA = 2 
    q10_Km_BicA = 2 
    q10_kin_Sco = 0.6 #Uemura paper

    # -- calculate conductivity of CO2 and bicarb through several concentric thylakoid layers -- 

    PlipCO2 = PlipCO2 * q10_PlipCO2**((temperature-25)/10) #um/s (CO2 permeability coefficient of a double lipid layer membrane)
    PlipHCO3 = PlipHCO3 * (1e6) * q10_PlipHCO3**((temperature-25)/10) #m/s, converted to um/s (bicarb permeability coefficient of a double lipid layer membrane)

    #generate list of surface areas of the stacked thylakoids (assume negligible thickness of the thylakoid membranes)

    thylakoid_areas = np.zeros(len(range(Nmem)))
    if Nmem != 1:
        thylakoid_additional_radius = (radCplast - radPyr) / (Nmem-1)
    else: 
        thylakoid_additional_radius = 0
    for i in range(Nmem):
        thylakoid_radius = radCplast - (i * thylakoid_additional_radius)
        thylakoid_areas[i] = 4*np.pi*(thylakoid_radius**2)
    
    #calculate total permeability of thylakoids TO CO2
    individual_resistance_thy_CO2 = np.zeros(len(thylakoid_areas))
    individual_resistance_total_CO2 = np.zeros(len(thylakoid_areas))
    for i in range(len(thylakoid_areas)):
        individual_resistance_thy_CO2[i] = 1/(PlipCO2*thylakoid_areas[i])
        individual_resistance_total_CO2[i] = 1/(PlipCO2*thylakoid_areas[i])
    individual_resistance_total_CO2[-1] = 1/(PlipCO2*SA_Cell)
    permeability_total_thy_CO2 = (sum(individual_resistance_thy_CO2))**(-1) #um^3/s
    permeability_total_cell_CO2 = (sum(individual_resistance_total_CO2))**(-1)

    #calculate total permeability of thylakoids TO HCO3-
    individual_resistance_thy_bicarb = np.zeros(len(thylakoid_areas))
    individual_resistance_total_bicarb = np.zeros(len(thylakoid_areas))
    for i in range(len(thylakoid_areas)):
        individual_resistance_thy_bicarb[i] = 1/(PlipHCO3*thylakoid_areas[i])
    permeability_total_thy_bicarb = (sum(individual_resistance_thy_bicarb))**(-1) #um^3/s

    #adjust oxygen permeability as described in Fridlyand paper
    PthyO2 = np.sqrt(44/32)*permeability_total_cell_CO2 #in the Fridlyand paper Po was 3.22e3 um^3/s

    # -- set all other parameters and adjust for temperature effects -- 

    #calculated with PlipCO2 or PlipHCO3, which are temperature-adjusted above (do not need further temperature adjustment)
    Pmc = (PlipCO2*SA_Cell) / 1e15            # um^3/s, converted to L/s (coefficient of CO2 permeability from medium to cytoplasm)
    Pcp = permeability_total_thy_CO2 / 1e15   # um^3/s, converted to L/s (coefficient of CO2 permeability from stroma to cytoplasm)
    Po = PthyO2 / 1e15                    # um^3/s, converted to L/s (coefficient of O2 permeability from thylakoid system to external medium)
    Pup = permeability_total_thy_bicarb / 1e15 # um^3/s, converted to L/s (coefficient of bicarb permeability from chloroplast inner border to stroma)
    Puc = (PlipHCO3*SA_Cell) / 1e15 # um^3/s, converted to L/s (coefficient of bicarb permeability from medium to cytoplasm)

    #other parameters with temperature adjustment already included or not relevant (do not need further temperature adjustment)
    CA_conc = CA_conc *1e6 / 1000  #mol m^-3, converted to uM, carbonic anhydrase concentration (using cytosolic value)
    pKa = 6            #CO2 to bicarb system overall pKa
    
    H_conc_cyt = 10**(-1*pH_cytosol)
    Keq_cyt = pH_cytosol - pKa 
    Keq_cyt = (10**Keq_cyt)*(H_conc_cyt)
    Keq_cyt = Keq_cyt * 1e6

    H_conc_cpl = 10**(-1*pH_chloroplast)
    Keq_cpl = pH_chloroplast - pKa
    Keq_cpl = (10**Keq_cpl)*(H_conc_cpl)
    Keq_cpl = Keq_cpl * 1e6

    protonc = H_conc_cyt*1e6  # uM
    protonp = H_conc_cpl*1e6  # uM 
    
    Henry_CO2 = 0.035 #mol kg^-1 bar^-1 (standard-temperature Henry's law constant, CO2) (use this rather than 0.034 to be consistent with 1st paper)
    Henry_temp_CO2 = 2400 #K (Henry's law temperature dependence constant, CO2)
    Henry_O2 = 0.0013 #mol kg^-1 bar^-1 (standard-temperature Henry's law constant, O2)
    Henry_temp_O2 = 1700 #K (Henry's law temperature dependence constant, O2)
    atm_pressure = 1.01325 #bar (atmospheric pressure)
    ppm_CO2 = 400/1000000 #ppm (the ppm of CO2 in the air)
    percent_O2 = 21 #% (the % of O2 in the air)

    #parameters that do need temperature adjustment (inorganic carbon concentrations)

    Henry_adjusted_CO2 = Henry_CO2*np.exp(Henry_temp_CO2*((1/(temperature+273.15))-(1/298.15))) #mol kg^-1 bar^-1 (temperature-adjusted Henry's law constant, CO2)
    Henry_adjusted_O2 = Henry_O2*np.exp(Henry_temp_O2*((1/(temperature+273.15))-(1/298.15))) #mol kg^-1 bar^-1 (temperature-adjusted Henry's law constant, O2)
    CO2m_molkg = Henry_adjusted_CO2 * (atm_pressure*ppm_CO2) #mol kg^-1 (medium [CO2])
    O2m_molkg = Henry_adjusted_O2 * (atm_pressure*(percent_O2/100)) #mol kg^-1 (medium [O2])
    water_density = ((-5e-6*(temperature**2)) + (8e-6*temperature) + 1.0001)*10e2*0.001 #kg m^-3, converted to kg L^-1 (density of water at temperature)

    CO2m_calc = CO2m_molkg*water_density*1e6 # uM (medium [CO2])
    O2m_calc = O2m_molkg*water_density*1e6  # uM (medium [O2])
    HCm_calc = 0              # uM (medium [bicarbonate])

    #other parameters that do need temperature adjustment
    Vc = Vc # mol/s (maximum rate of CO2 fixation inside the chloroplast)
    Kc = Kc  # uM (apparent Michaelis constant for CO2 binding in the Calvin/CBB cycle)
    Ko = Ko
    k2 = 6.2e-2 * q10_k2**((temperature-25)/10)   #s^-1
    k_2 = 1.37e1 * q10_k_2**((temperature-25)/10)  #s^-1
    ShCO2_forward = 6e-2 * q10_ShCO2_forward**((temperature-25)/10) # s^-1 # spontaneous hydration of CO2 (ShCO2)
    ShCO2_reverse = 2e1 * q10_ShCO2_reverse**((temperature-25)/10) # s^-1
    khcat = khcat * q10_khcat**((temperature-25)/10) #s^-1, hydration rate constant
    KCO2 = KCO2 * 1e6 /1000 * q10_KCO2**((temperature-25)/10) #mol m^-3, converted to uM, affinity of CA for CO2
    KHC = KHC * 1e6 /1000 * q10_KHC**((temperature-25)/10) #mol m^-3, converted to uM, affinity of CA for bicarbonate
    Vmax_BicA = ((Vmax_BicA)/1e12) * q10_Vmax_BicA**((temperature-25)/10) # mol um^-2 s^-1, converted to...?
    Km_BicA = Km_BicA * q10_Km_BicA**((temperature-25)/10) # umol L^-1
    kin_Sco = kin_Sco * q10_kin_Sco**((temperature-25)/10)     #a ratio (Uemura paper)
    
    DiffusionCoefficientO2 = 3.05e-5 * 1e8
    DiffusionCoefficientCO2 = (DiffusionCoefficientO2 / (np.sqrt(32/44))) #m^2 s

    Pmb = ((SA_Boundary * DiffusionCoefficientCO2) / (LengthBoundary)) / 1e15
    Pmo = ((SA_Boundary * DiffusionCoefficientO2) / (LengthBoundary)) / 1e15
       
    return T1, T2, T3, Pmc, Pcp, Pup, Puc, Po, CO2m_calc, O2m_calc, HCm_calc, CA_conc, Vc, Kc, Ko, k2, k_2, ShCO2_forward, ShCO2_reverse, khcat, KCO2, KHC, Vmax_BicA, Km_BicA, kin_Sco, Keq_cyt, Keq_cpl, Pmb, T4, pH_chloroplast, pH_cytosol, q10_PlipCO2, q10_PlipHCO3, PlipHCO3, pumpCost, PlipCO2, radCell, protonc, protonp, SA_Cplast, Nmem, Pmo

filename_number = sys.argv[-1]
filename = 'Testparamssubset_{0}.csv'.format(sys.argv[-1])
input_dataframe = pd.read_csv(filename,index_col=0)


data = {
    'T1':[],
    'T2':[],
    'T3':[],
    'Pmc':[],
    'Pcp':[],
    'Pup':[],
    'Puc':[],
    'Po':[],
    'CO2m':[],
    'O2m':[],
    'HCm':[],
    'CA_conc':[],
    'Vc':[],
    'Kc':[],
    'Ko':[],
    'k2':[],
    'k_2':[],
    'ShCO2_forward':[],
    'ShCO2_reverse':[],
    'khcat':[],
    'KCO2':[],
    'KHC':[],
    'Vmax_BicA':[],
    'Km_BicA':[],
    'kin_Sco':[],
    'Keq_cyt':[], 
    'Keq_cpl':[],
    'Pmb':[],
    'T4':[],
    'pH_chloroplast':[],
    'pH_cytosol':[],
    'q10_Plip_CO2':[],
    'q10_Plip_HCO3':[],
    'PlipHCO3':[],
    'pumpCost':[],
    'PlipCO2':[],
    'radCell':[],
    'protonc':[],
    'protonp':[],
    'SA_Cplast':[],
    'Stacks':[],
    'Pmo':[],
    'Compensation_point':[],
    'CO2 cytosol':[],
    'CO2 boundary:':[],
    'HCO3 cytosol:':[],
    'CO2 chloroplast':[],
    'HCO3 chloroplast':[],
    'O2 chloroplast':[],
    'O2 boundary':[],
    'Net Assimilation':[],
    'VO/VC':[],
    'ATP value':[],
    'NADPH value':[],
    'Leak %':[],
    'ATP per CO2':[],
    'Steady State?':[]}


my_dataframe = pd.DataFrame(data)

### Add in loops for each additional parameter you're testing ###

for x in range(len(input_dataframe)):
    arguments = input_dataframe

    root, sol_1, sol_2, sol_3, sol_4, sol_5, sol_6, sol_7, net_assim, vo_vc, ATP_value, NADPH_value, leak, ATP_per_CO2, steady_state = compensation_calculation(arguments.iloc[x])

    my_dataframe.loc[len(my_dataframe.index)] = [arguments.iloc[x][0], arguments.iloc[x][1], arguments.iloc[x][2], arguments.iloc[x][3], arguments.iloc[x][4],
                                                    arguments.iloc[x][5], arguments.iloc[x][6], arguments.iloc[x][7], arguments.iloc[x][8], arguments.iloc[x][9],
                                                    arguments.iloc[x][10], arguments.iloc[x][11], arguments.iloc[x][12], arguments.iloc[x][13], arguments.iloc[x][14],
                                                    arguments.iloc[x][15], arguments.iloc[x][16], arguments.iloc[x][17], arguments.iloc[x][18], arguments.iloc[x][19],
                                                    arguments.iloc[x][20], arguments.iloc[x][21], arguments.iloc[x][22], arguments.iloc[x][23], arguments.iloc[x][24],
                                                    arguments.iloc[x][25], arguments.iloc[x][26], arguments.iloc[x][27], arguments.iloc[x][28], arguments.iloc[x][29], 
                                                    arguments.iloc[x][30], arguments.iloc[x][31], arguments.iloc[x][32], arguments.iloc[x][33], arguments.iloc[x][34], 
                                                    arguments.iloc[x][35], arguments.iloc[x][36], arguments.iloc[x][37], arguments.iloc[x][38], arguments.iloc[x][39],arguments.iloc[x][40],arguments.iloc[x][41], root,
                                                    sol_1, sol_2, sol_3, sol_4, sol_5, sol_6, sol_7, net_assim, vo_vc, ATP_value, NADPH_value, leak, ATP_per_CO2, steady_state]
    
output_filename = 'Results_{0}.csv'.format(sys.argv[-1])    

my_dataframe.to_csv(output_filename)   