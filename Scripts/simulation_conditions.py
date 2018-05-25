# -*- coding: utf-8 -*-
"""
Creates files with the simulation conditions
"""

from __future__ import division
import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import minimize
import argparse

#Before running script run, "pip install pymbar, pip install CoolProp"

REFPROP_path='/home/ram9/REFPROP-cmake/build/' #Change this for a different system

CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH,REFPROP_path)

### Physical constants
N_A = 6.0221409e23 #Avogadro's number [1/mol]
conv_m3tonm3 = 1e-27 #[m3/nm3]

compound_dic = {'Propane':'C3H8','Butane':'C4H10','Octane':'C8H18','Isobutane':'IC4H10','Isopentane':'IC5H12','Isohexane':'IC6H14','Isooctane':'IC8H18','Neopentane':'NEOC5H12'}

def convert_rhol_Lbox(rhol,Nmol,Mw):
    nrho = rhol/Mw*N_A*conv_m3tonm3 #[1/nm3]
    Vbox= Nmol / nrho #[nm3]
    Lbox = Vbox**(1./3.)
    return Lbox

def calc_Tsat_rhol(rhol,TPT,TC,compound):
    Tsat_guess = np.average([TPT,TC])
    dev = lambda Tsat: (CP.PropsSI('D','T',Tsat,'Q',0,'REFPROP::'+compound) - rhol)**2.
    opt = minimize(dev,Tsat_guess)
    return np.round(opt.x)

def extrapolate_rho(press,rho,press_range,rho_range):
    press_valid = press_range[np.isfinite(rho_range)]
    rho_valid = rho_range[np.isfinite(rho_range)]

    press_high = press_valid[-int(len(press_valid)/5):]
    rho_high = rho_valid[-int(len(rho_valid)/5):]

    pfit = np.polyfit(press_high,rho_high,deg=1)

    rho[np.isinf(rho)] = np.polyval(pfit,press[np.isinf(rho)]) #Only extrapolate for the densities that are infinite

    return rho     

def create_files(Temp_sim,rhol_sim,Lbox_sim,Nmol,filepath,press_sim=[]):

    f = open(filepath+'Nmol','w')
    f.write(str(Nmol))
    f.close()
    
    f = open(filepath+'Temp','w')
    for Temp in Temp_sim:
        f.write(str(Temp)+' \n')
    f.write('Endline')
    f.close()
    
    f = open(filepath+'rho','w')
    for rho in rhol_sim:
        f.write(str(rho)+' \n')
    f.write('Endline')
    f.close()
    
    f = open(filepath+'liquid_box_N'+str(Nmol),'w')
    for Lbox in Lbox_sim:
        f.write(str(Lbox)+' \n')
    f.write('Endline')
    f.close()
    
    try:
        if len(press_sim) > 0:
            f = open(filepath+'press','w')
            for press in press_sim:
                f.write(str(press)+' \n')
            f.write('Endline')
            f.close()
    except:
        pass

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--comp",type=str,help="Specify the compound to analyze")
    parser.add_argument("-N","--Nmol",type=int,help="Specify the number of moles")
    parser.add_argument("--sat",action='store_true',help="If saturation conditions are desired")
    parser.add_argument("--T293highP",action='store_true',help="If T=293 K and high pressure conditions are desired")
    parser.add_argument("--nstates",type=int,help="Specify the number of state points")
    
    args = parser.parse_args()
    
    if args.nstates:
        nstates = args.nstates
    else:
        nstates = 5
    
    if args.comp and args.Nmol:

        compound=args.comp
        Nmol = args.Nmol
        
        Mw = CP.PropsSI('M','REFPROP::'+compound) #[kg/mol]
        RP_TC = CP.PropsSI('TCRIT','REFPROP::'+compound) #[K]
        RP_Tmin =  CP.PropsSI('TMIN','REFPROP::'+compound) #[K]
        RP_Pmax = CP.PropsSI('PMAX','REFPROP::'+compound)/(10**6) #[MPa]
          
        if args.sat:
            
            if nstates > 5:
                Tmin = RP_Tmin
            else:
                Tmin = 0.45 * RP_TC

            Tmax = 0.85 * RP_TC

            Tsat_sim = np.trunc(np.linspace(Tmin,Tmax,num=nstates))
            rhol_sim = CP.PropsSI('D','T',Tsat_sim,'Q',0,'REFPROP::'+compound) #[kg/m3]
            Lbox_sim = convert_rhol_Lbox(rhol_sim,Nmol,Mw)
            
            file_path = 'Saturation_Conditions/'+compound_dic[compound]+'_'
            
            create_files(Tsat_sim,rhol_sim,Lbox_sim,Nmol,file_path)
        
        if args.T293highP:
            
            Tsat_sim = np.array([293.]*nstates)

            Psat = CP.PropsSI('P','T',Tsat_sim[0],'Q',0,'REFPROP::'+compound)/ (10**6) #[MPa]        

            ### Although REFPROP provide Pmax (RP_Pmax above), this value appears to depend on the temperature.
            ### For this reason, we determine the maximum pressure independently

            press_range = np.trunc(np.linspace(0,1000,10000)) #[MPa]
            press_range[0] = np.max([0.1,Psat]) #First value should be either 0.1 MPa or Psat so that it is liquid
            press_range_Pa = press_range*10**6 #[Pa], REFPROP requires Pa
            
            rho_range = CP.PropsSI('D','T',Tsat_sim[0],'P',press_range_Pa,'REFPROP::'+compound) #[kg/m3]
            visc_range = CP.PropsSI('VISCOSITY','T',Tsat_sim[0],'P',press_range_Pa,'REFPROP::'+compound)

            press_max = np.max(press_range[np.isfinite(rho_range)]) #We do not want the points above where REFPROP provides density
            press_max_visc = np.max(press_range[np.isfinite(visc_range)]) #The pressure max for viscosity is often different (although this depends on temperature)
            #print(press_max,press_max_visc,RP_Pmax)

            ### With the maximum pressure we repeat this process but only with nstates
            ### We must use press_max, not press_max_visc in this case

#            press_sim = np.trunc(np.linspace(0,press_max,nstates)) #[MPa]
#            press_sim[0] = np.max([0.1,Psat]) #First value should be either 0.1 MPa or Psat so that it is liquid
#            press_sim_Pa = press_sim*10**6 #[Pa], REFPROP requires Pa

#            rhol_sim = CP.PropsSI('D','T',Tsat_sim[0],'P',press_sim_Pa,'REFPROP::'+compound) #[kg/m3]

            ### Alternatively we can try extrapolating to high pressure
            ### Here we can use press_max_visc because we will just extrapolate rho in this case

            press_sim = np.trunc(np.linspace(0,1000,nstates)) #[MPa]
            press_sim[0] = np.max([0.1,Psat]) #First value should be either 0.1 MPa or Psat so that it is liquid
            press_sim[1] = np.min([press_sim[1],press_max_visc])

            press_sim_Pa = press_sim*10**6 #[Pa], REFPROP requires Pa

            rhol_sim = CP.PropsSI('D','T',Tsat_sim[0],'P',press_sim_Pa,'REFPROP::'+compound) #[kg/m3]

            rhol_sim = extrapolate_rho(press_sim,rhol_sim,press_range,rho_range)

            press_sim = press_sim*10. #[bar] for Gromacs
            Lbox_sim = convert_rhol_Lbox(rhol_sim,Nmol,Mw)
            
            file_path = 'T293highP_Conditions/'+compound_dic[compound]+'_'
            
            create_files(Tsat_sim,rhol_sim,Lbox_sim,Nmol,file_path,press_sim)

if __name__ == '__main__':
    '''
    python simulation_conditions.py --comp str --Nmol int --sat --T293highP
  
    '''

    main()  