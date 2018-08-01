"""
This script is used to find data in a file hierarchy described by
/compound/Gromacs/[Condition]_Viscosity/[Forcefield]_N[Number]_LINCS
This can fairly easily be changed.
Files it loads are SaturationSettings/SaturatedSettings or T293highPSettings
Lbox_all, press_all_log, GK_eta_boots_rho# and GK_MCMC_avg_rho#
"""

from __future__ import division
import numpy as np 
import argparse
import os
import math

N_A = 6.022e23
nm3tom3 = 1e-27
gmtokg = 1e-3

# This function was obtained from randlet.com/blog/python-significant-figures-format/
# which is also available on github.
# It prints a number to the precision desired.
def to_precision(x,p):
    """
    returns a string representation of x formatted with a precision of p

    Based on the webkit javascript implementation taken from here:
    https://code.google.com/p/webkit-mirror/source/browse/JavaScriptCore/kjs/number_object.cpp
    """

    x = float(x)

    if x == 0.:
        return "0." + "0"*(p-1)

    out = []

    if x < 0:
        out.append("-")
        x = -x

    e = int(math.log10(x))
    tens = math.pow(10, e - p + 1)
    n = math.floor(x/tens)

    if n < math.pow(10, p - 1):
        e = e -1
        tens = math.pow(10, e - p+1)
        n = math.floor(x / tens)

    if abs((n + 1.) * tens - x) <= abs(n * tens -x):
        n = n + 1

    if n >= math.pow(10,p):
        n = n / 10.
        e = e + 1

    m = "%.*g" % (p, n)

    if e < -2 or e >= p:
        out.append(m[0])
        if p > 1:
            out.append(".")
            out.extend(m[1:p])
        out.append('e')
        if e > 0:
            out.append("+")
        out.append(str(e))
    elif e == (p -1):
        out.append(m)
    elif e >= 0:
        out.append(m[:e+1])
        if e+1 < len(m):
            out.append(".")
            out.extend(m[e+1:])
    else:
        out.append("0.")
        out.extend(["0"]*-(e+1))
        out.append(m)

    return "".join(out)

# Determine the number of significant figures to print
# given a number and its uncertainty
def det_digits(number, uncer):
    # Who cares about 0 anyway?
    if number != 0:
        number = math.log10(abs(number))
    if uncer != 0:
        uncer = math.log10(abs(uncer))
    number = math.floor(number)
    uncer = math.floor(uncer)
    return int(number - uncer + 2)

# Convenience wrapper for det_digits and to_precision to 
# convert to the proper precision in one step based on
# uncertainties
def sig_figs(number, uncer):
    digits = det_digits(number, uncer)
    return str(to_precision(number, digits))

def Lbox_to_rho(Lbox,Nsim,Mw):
    '''
     Lbox in nm, Mw in gm/mol
     Returns rho in kg/m3
    '''
    Vol = Lbox**3. #[nm3/N]
    Vol /= Nsim #[nm3/molecule]
    Vol *= N_A #[nm3/mol]
    Vol /= Mw #[nm3/gm]
    Vol *= nm3tom3 #[m3/gm]
    Vol *= 1000. #[m3/kg]
    rho = 1./Vol #[kg/m3]
    return rho 

MW_list = {'C2H6':30.07,'C3H8':44.096,'C4H10':58.122,'C8H18':114.23,'C12H26':170.33,'IC4H10':58.122,'IC8H18':114.23,'IC5H12':72.149,'3MPentane':86.2,'23DMButane':86.1754,'NEOC5H12':72.149,'C16H34':226.41,'C22H46':310.6027,'IC6H14':86.175}

# Read in the bootstrap and average values from a given filepath and rho
# if they exist
# Returns an empy array if they do not
def read_boots(filepath, rho):
    try:
        bootsdata = np.loadtxt(filepath+"GK_eta_boots_rho"+str(rho))
    except:
        return []
    return bootsdata

# Takes a bootstrap list and returns the mean, low 95% bound, and high 95% bond
def process_boots(bootsdata):
    bootsdata = np.sort(bootsdata)
    i95 = int(0.025*len(bootsdata))
    bootsdata_95 = bootsdata[i95:-i95]
    avg = np.mean(bootsdata_95)
    low = avg - bootsdata_95[0]
    high = bootsdata_95[-1] - avg
    mean = np.mean(bootsdata)
    return mean, mean - low, mean + high 

# Load all of the data for a given force field, compound, and condition
# rhos are index values, indicating which boots files were actually loaded
# when a boots file is skipped, that rho is not added to the rhos list and
# a place holder "0" is added to the other arrays to avoid index errors when
# data are printed in coordination with the settings files' information
def load_forcefield_data(compound, condition, model, number):
    root_path = compound+"/Gromacs/"+condition+"_Viscosity/"+model+"_N"+str(number)+"_LINCS/"
    done = False
    rho = 0
    means = []
    lows = []
    highs = []
    # Because a rho may be absent, this needs to be kept track of
    rhos = []
    # Until a load fails, continue trying to load the next rho
    # If some rho under 5 does not exist, keep trying to load data
    while done == False or rho < 5:
        boots = read_boots(root_path,rho)
        if len(boots) > 0:  # Load succeeded
            mean, low, high = process_boots(boots)
            means.append(mean)
            lows.append(low)
            highs.append(high)
            rhos.append(rho)
        else:  # Load failed
            means.append(0)
            lows.append(0)
            highs.append(0)
            rhos.append(rho)
            done = True
        rho = rho + 1
    # An extra rho/mean/low/high has been appended here, we take it off
    means = means[:-1]
    lows = lows[:-1]
    highs = highs[:-1]
    rhos = rhos[:-1]
    return means, lows, highs, rhos

# Load the simulation conditions for a given compound, condition, and field
def load_sim_conditions_sat(compound, model, number):
    root_path = compound+"/Gromacs/Saturation_Viscosity/"+model+"_N"+str(number)+"_LINCS/"
    # Load code taken cfrom compare_force_fields.py
    try:
        #print("trying "+root_path+'SaturatedSettings.txt')
        sim_sat = np.loadtxt(root_path+'SaturatedSettings.txt',skiprows=1)
        #print("Found it!")
    except:
        sim_sat = np.loadtxt(root_path+'SaturationSettings.txt',skiprows=1)
    try:
        Tsat_sim = sim_sat[:,2]
        L_sat_sim = sim_sat[:,1]
    except:
        Tsat_sim = np.array([sim_sat[2]])
        L_sat_sim = np.array([sim_sat[1]])

    density = Lbox_to_rho(L_sat_sim,number,MW_list[compound])
    return Tsat_sim, density

# Returns five arrays of data and uncertainties because I am
# too lazy to find out how to concatenate them into one array
# Returns: temperature, liquid box size, pressure, box uncertainty,
# pressure uncertainty
# Returns an empty value and prints a warning on an error
def load_sim_conditions_highP(compound, model, number):
    root_path = compound+"/Gromacs/T293highP_Viscosity/"+model+"_N"+str(number)+"_LINCS/"
    try:
        data = np.loadtxt(root_path+"T293highPSettings.txt",skiprows=1)
    except Exception, e:
        print("Warning: failed to load data from "+root_path+": "+str(e))
        return [], [], [], [] 
    t_sim = []
    press_sim = []
    upress_sim = []
    u_density = []
    
    try:
        t_sim = data[:,2]
        L_sim = data[:,1]
    except:
        t_sim = np.array([data[2]])
        L_sim = np.array([data[1]])
    # This variable will be reused if we load an Lbox file
    density = Lbox_to_rho(L_sim,number,MW_list[compound])

    nrho = len(t_sim)
    nrhomax = nrho  # Currently redundant; may be an argument later
    try:
        press_sim = np.loadtxt(root_path+"press_all_log")  # Bar
        if len(press_sim) != 60*nrho:
            print("Warning on "+root_path+": length of pressure log != 60*nrho")
        press_sim /= 10.
        press_sim = press_sim[:60*nrho]
        # Change it into a number of arrays evenly distributing the pressures
        press_sim = press_sim.reshape([nrho,int(len(press_sim)/nrho)])
        upress_sim = 1.96*np.std(press_sim,axis=1)  # Take stdev
        #print(compound+" upress_sim")
        #print(upress_sim)
        #print("done")
        press_sim = np.mean(press_sim,axis=1)  # compute mean pressures
    except Exception, e:
        print("No press_all_log file on "+root_path+": "+str(e))
    # Try to load Lbox to create a density vector
    try:
        L_sim = np.loadtxt(root_path+'Lbox_all') #[bar]
        # Attempt to compute a density vector from the Lbox values
        density = Lbox_to_rho(L_sim,number,MW_list[compound])
        if len(density) != 60*nrho:
            print("Warning on "+root_path+": length of Lbox log != 60*nrho (nrho ="+str(nrho)+")")
        density = density[:60*nrhomax]
        density = density.reshape([nrho,int(len(density)/nrho)])
        u_density = 1.96*np.std(density,axis=1)
        density = np.mean(density,axis=1)
    except Exception, e:
        print("No Lbox file on "+root_path+": "+str(e))

    return t_sim, density, press_sim, upress_sim, u_density

# Print loaded saturation information to a file
def print_sat_info(compound, model, tsat, dens_sim, satmeans, satlows, sathighs, satrhos):
    # For each rho which was actually available
    # Apply header
    sf = open("compiled_data/"+compound+"_sat_"+model+".txt","w")
    sf.write("T (K) 	rho (Kg/m^3)	eta (Pa-s)	low eta 95% (Pa-s)	high eta 95% (Pa-s)\n")
    for rho in satrhos:  
        sf.write(str(tsat[rho])+"	"+str(dens_sim[rho])+"	"+str(satmeans[rho])+"	"+str(satlows[rho])+"	"+str(sathighs[rho]))
        sf.write("\n")
    sf.close()

# Print high pressure information into a file in one of two formats depending on pressure data availability
def print_T293_info(compound, model, t_sim, dens_sim, press_sim, upress_sim, hpmeans, hplows, hphighs, hprhos):
    hf = open("compiled_data/"+compound+"_highP_"+model+".txt","w")
   
    # Apply appropriate header depending on presence of pressure information
    if len(press_sim) > 0:
        hf.write("T (K)	P (MPa)	Uncertainty (MPa)	eta (Pa-s)	low eta 95% (Pa-s)	high eta 95% (Pa-s)\n")
    else:
        hf.write("T (K)	rho (Kg/m^3)	eta (Pa-s)	low eta 95% (Pa-s)	high eta 95% (Pa-s)\n")

    for rho in hprhos:
        if len(press_sim) > 0:  # pressure values avaialble for T293
            hf.write(str(t_sim[rho])+"	"+str(press_sim[rho])+"	"+str(upress_sim[rho])+"	"+str(hpmeans[rho])+"	"+str(hplows[rho])+"	"+str(hphighs[rho]))
        else:  # No pressure values available
            hf.write(str(t_sim[rho])+"	"+str(dens_sim[rho])+"	"+str(hpmeans[rho])+"	"+str(hplows[rho])+"	"+str(hphighs[rho]))
        hf.write("\n")
    hf.close()

# Iterate through means, lows, and highs and determine the
# maximum uncertainty represented for each point.
def determine_uncertainties(means, lows, highs):
    uncers = []
    for i in range(len(means)):
        maximum = max(means[i] - lows[i], highs[i] - means[i])
        uncers.append(maximum)
    return uncers

def handle_T293_info(comp, number):
    # Flags are used to determine overall success
    failed_P = False
    failed_TA = False
    failed_Tr = False
    # These are used to determine how many state points need to be printed
    len_P = 0
    len_TA = 0
    len_Tr = 0

    # As these lengths must be referenced later, they are declared here
    TA_udens = []
    P_udens = []
    Tr_udens = []

     
    # Load forcefield data and determine maximum uncertainties for each of the
    # Compounds in turn. Following that, open a file and write the information
    # hprhos are a relic and can safely be ignored
    try:
        P_hpmeans, P_hplows, P_hphighs, P_hprhos = load_forcefield_data(comp,"T293highP","Potoff",number)
        P_t_sim, P_dens_sim, P_press_sim, P_upress_sim, P_udens = load_sim_conditions_highP(comp,"Potoff",number)
        P_uncer = determine_uncertainties(P_hpmeans, P_hplows, P_hphighs)
        len_P = len(P_hpmeans)
    except Exception, e:
        failed_P = True
        print("Unable to load T293 Potoff info for "+comp+": "+str(e))
    try:
        TA_hpmeans, TA_hplows, TA_hphighs, TA_hprhos = load_forcefield_data(comp,"T293highP","TAMie",number)
        TA_t_sim, TA_dens_sim, TA_press_sim, TA_upress_sim, TA_udens = load_sim_conditions_highP(comp,"TAMie",number)
        TA_uncer = determine_uncertainties(TA_hpmeans, TA_hplows, TA_hphighs)
        len_TA= len(TA_hpmeans)
    except Exception, e:
        failed_TA = True
        print("Unable to load T293 TAMie info for "+comp+": "+str(e))
    try:
        Tr_hpmeans, Tr_hplows, Tr_hphighs, Tr_hprhos = load_forcefield_data(comp,"T293highP","TraPPE",number)
        Tr_t_sim, Tr_dens_sim, Tr_press_sim, Tr_upress_sim, Tr_udens = load_sim_conditions_highP(comp,"TraPPE",number)
        Tr_uncer = determine_uncertainties(Tr_hpmeans, Tr_hplows, Tr_hphighs)
        len_Tr = len(Tr_hpmeans)
    except Exception, e:
        failed_Tr = True
        print("Unable to load T293 TraPPE info for "+comp+": "+str(e))

    points = 0  # Number of points to print eventually
    if failed_P and failed_TA and failed_Tr:
    # No data could be loaded
        return
    else:
        points = max(len_P, len_TA, len_Tr)    

    hf = open("compiled_data/"+comp+"_T293highP.txt","w")
        

    # If Lbox informaiton/density uncertainty info is available
    if len(P_udens) > 0 or len(TA_udens) > 0 or len(Tr_udens) > 0:
        # sort data into the proper pressure order if possible
        try:
            if not failed_P:
                P_press_sim, P_dens_sim, P_udens, P_upress_sim, P_hpmeans, P_uncer  = zip(*sorted(zip(P_press_sim,P_dens_sim,P_udens,P_upress_sim,P_hpmeans,P_uncer))) 
            if not failed_TA:
                TA_press_sim, TA_dens_sim, TA_udens, TA_upress_sim, TA_hpmeans, TA_uncer  = zip(*sorted(zip(TA_press_sim,TA_dens_sim,TA_udens,TA_upress_sim,TA_hpmeans,TA_uncer))) 
            if not failed_Tr:
                Tr_press_sim, Tr_dens_sim, Tr_udens, Tr_upress_sim, Tr_hpmeans, Tr_uncer  = zip(*sorted(zip(Tr_press_sim,Tr_dens_sim,Tr_udens,Tr_upress_sim,Tr_hpmeans,Tr_uncer))) 
        except Exception, e:
            print("Unable to sort all T293 data for "+comp+": "+str(e))
       
        hf.write("	Potoff			TAMie			TraPPE\n")
        hf.write("$T$ [K]	$\\rho^{\\rm comp}_{\\rm liq}$ [Kg/m$^3$]	$P$ [MPa]	$\eta^{\\rm comp}_{\\rm liq}$ [Pa-s]	$\\rho^{\\rm comp}_{\\rm liq}$ [Kg/m$^3$]	$P$ [MPa]	$\eta^{\\rm comp}_{\\rm liq}$ [Pa-s]	$\\rho^{\\rm comp}_{\\rm liq}$ [Kg/m$^3$]	$P$ [MPa]	$\eta^{\\rm comp}_{\\rm liq}$ [Pa-s]\n")
        # Try to sort data into the proper pressure order
  
        for rho in range(points):
            hf.write("{0:.3g}	".format(293))
            # The two checks here are just becaue they usually catch things and I don't want to have to write checks for every array
            if failed_P or rho >= len(P_hprhos) or rho >= len(P_udens):
                hf.write("X	X	X	")
            else:
                hf.write(sig_figs(P_dens_sim[rho],P_udens[rho])+"$_{"+to_precision(P_udens[rho],2)+"}$	"+sig_figs(P_press_sim[rho],P_upress_sim[rho])+"$_{"+to_precision(P_upress_sim[rho],2)+"}$	"+sig_figs(P_hpmeans[rho],P_uncer[rho])+"$_{"+to_precision(P_uncer[rho],2)+"}$	")
            if failed_TA or rho >= len(TA_hprhos) or rho >= len(TA_udens):  # We've run out of room or skipped this one
                hf.write("X	X	X	")
            else:
                hf.write(sig_figs(TA_dens_sim[rho],TA_udens[rho])+"$_{"+to_precision(TA_udens[rho],2)+"}$	"+sig_figs(TA_press_sim[rho],TA_upress_sim[rho])+"$_{"+to_precision(TA_upress_sim[rho],2)+"}$	"+sig_figs(TA_hpmeans[rho],TA_uncer[rho])+"$_{"+to_precision(TA_uncer[rho],2)+"}$	")
            if failed_Tr or rho >= len(Tr_hprhos) or rho >= len(Tr_udens):  # We've run out of room
                hf.write("X	X	X\n")
            else:
                hf.write(sig_figs(Tr_dens_sim[rho],Tr_udens[rho])+"$_{"+to_precision(Tr_udens[rho],2)+"}$	"+sig_figs(Tr_press_sim[rho],Tr_upress_sim[rho])+"$_{"+to_precision(Tr_upress_sim[rho],2)+"}$	"+sig_figs(Tr_hpmeans[rho],Tr_uncer[rho])+"$_{"+to_precision(Tr_uncer[rho],2)+"}$\n")

    # No Lbox data available; use default density
    else:
        try:
            if not failed_P:
                P_press_sim, P_dens_sim, P_hpmeans, P_uncer = zip(*sorted(zip(P_press_sim,P_dens_sim,P_hpmeans,P_uncer))) 
            if not failed_TA:
                TA_press_sim, TA_dens_sim, TA_hpmeans, TA_uncer = zip(*sorted(zip(TA_press_sim,TA_dens_sim,TA_hpmeans,TA_uncer))) 
            if not failed_Tr:
                Tr_press_sim, Tr_dens_sim, Tr_hpmeans, Tr_uncer = zip(*sorted(zip(Tr_press_sim,Tr_dens_sim,Tr_hpmeans,Tr_uncer))) 
        except Exception, e:
            print("Unable to sort T293 data for "+comp+": "+str(e))

        # Get some valid data into the P_dens_sim vector
        # and make sure the longest vector is employed!
        if failed_P or len_P < len_Tr or len_P < len_TA:
           if len_Tr > len_TA:
               P_dens_sim = Tr_dens_sim
           else:
               P_dens_sim = TA_dens_sim
        hf.write("		Potoff		TAMie		TraPPE\n")
        print("No Lbox data for T293highP "+comp)
         
        hf.write("$T$ [K]	$\\rho^{\\rm comp}_{\\rm liq}$ [Kg/m$^3$]	$P$ [MPa]	$\eta^{\\rm comp}_{\\rm liq}$ [Pa-s]	$P$ [MPa]	$\eta^{\\rm comp}_{\\rm liq}$ [Pa-s]	$P$ [MPa]	$\eta^{\\rm comp}_{\\rm liq}$ [Pa-s]\n")
        for rho in range(points):
            hf.write(to_precision(293,3)+"	"+to_precision(P_dens_sim[rho],5)+"	")
            if failed_P or rho >= len(P_hprhos):
                hf.write("X	X	")
            else:
                hf.write(sig_figs(P_press_sim[rho],P_upress_sim[rho])+"$_{"+to_precision(P_upress_sim[rho],2)+"}$	"+sig_figs(P_hpmeans[rho],P_uncer[rho])+"$_{"+to_precision(P_uncer[rho],2)+"}$	")
            if failed_TA or rho >= len(TA_dens_sim):
                hf.write("X	X	")
            else:
                hf.write(sig_figs(TA_press_sim[rho],TA_upress_sim[rho])+"$_{"+to_precision(TA_upress_sim[rho],2)+"}$	"+sig_figs(TA_hpmeans[rho],TA_uncer[rho])+"$_{"+to_precision(TA_uncer[rho],2)+"}$	")
            if failed_Tr or rho >= len(Tr_dens_sim):
                hf.write("X	X\n")
            else:
                hf.write(sig_figs(Tr_press_sim[rho],Tr_upress_sim[rho])+"$_{"+to_precision(Tr_upress_sim[rho],2)+"}$	"+sig_figs(Tr_hpmeans[rho],Tr_uncer[rho])+"$_{"+to_precision(Tr_uncer[rho],2)+"}$\n")


def handle_sat_info(comp, number):
    # Load forcefield data and determine the maximum uncertainties for each.
    # Flags are used to determine overall success
    failed_P = False
    failed_TA = False
    failed_Tr = False
    # These are used to determine how many state points need to be printed
    len_P = 0
    len_TA = 0
    len_Tr = 0
    # satrhos are a relic and can safely be ingored or used to pass in other information 
    try:
        P_t_sat, P_dens_sat = load_sim_conditions_sat(comp,"Potoff",number)
        P_satmeans, P_satlows, P_sathighs, P_satrhos = load_forcefield_data(comp,"Saturation","Potoff",number)
        P_uncer = determine_uncertainties(P_satmeans,P_satlows,P_sathighs)
        len_P = len(P_satmeans)
    except Exception, e:
        failed_P = True
        print("Cannot load saturation potoff data for "+comp+": "+str(e))
    try:
        TA_t_sat, TA_dens_sat = load_sim_conditions_sat(comp,"TAMie",number)
        TA_satmeans, TA_satlows, TA_sathighs, TA_satrhos = load_forcefield_data(comp,"Saturation","TAMie",number)
        TA_uncer = determine_uncertainties(TA_satmeans,TA_satlows,TA_sathighs)
        len_TA = len(TA_satmeans)
    except Exception, e:
        failed_TA = True
        print("Cannot load saturation TAMie data for "+comp+": "+str(e))
    try:
        Tr_t_sat, Tr_dens_sat = load_sim_conditions_sat(comp,"TraPPE",number)
        Tr_satmeans, Tr_satlows, Tr_sathighs, Tr_satrhos = load_forcefield_data(comp,"Saturation","TraPPE",number)
        Tr_uncer = determine_uncertainties(Tr_satmeans,Tr_satlows,Tr_sathighs)
        len_Tr = len(Tr_satmeans)
    except Exception, e:
        failed_Tr = True
        print("Cannot load saturation TraPPE data for "+comp+": "+str(e))
    try:
        if not failed_TA:
            TA_t_sat, TA_dens_sat, TA_satmeans, TA_uncer, TA_satrhos = zip(*sorted(zip(TA_t_sat, TA_dens_sat, TA_satmeans, TA_uncer, TA_satrhos)))
        if not failed_Tr:
            Tr_t_sat, Tr_dens_sat, Tr_satmeans, Tr_uncer, Tr_satrhos = zip(*sorted(zip(Tr_t_sat, Tr_dens_sat, Tr_satmeans, Tr_uncer, Tr_satrhos)))
        if not failed_P:
            P_t_sat, P_dens_sat, P_satmeans, P_uncer, P_satrhos = zip(*sorted(zip(P_t_sat, P_dens_sat, P_satmeans, P_uncer, P_satrhos)))
    except Exception, e:
        print("Sorting not possible for saturation "+comp+": "+str(e))

    # Dead; nothing loaded
    points = 0
    if failed_P and failed_TA and failed_Tr:
        return
    else: 
        points = max(len_P, len_TA, len_Tr)

    if failed_P or len_P < len_TA or len_P < len_Tr:
       # P_t_sat and P_dens_sat are printed into the file so it needs to exist
       # Make sure to get the longest vector to print
       if len_TA > len_Tr:
           P_t_sat = TA_t_sat
           P_dens_sat = TA_dens_sat
       else:
           P_t_sat = Tr_t_sat
           P_dens_sat = Tr_dens_sat

    # Open a file, loop over teh saturaiton points and write into the file
    sf = open("compiled_data/"+comp+"_sat.txt","w")
    sf.write("		Potoff	TAMie	TraPPE\n")
    sf.write("$T^{\\rm sat}$ [K]	$\\rho^{\\rm sat}_{\\rm liq}$ [Kg/m$^3$]	$\eta^{\\rm sat}_{\\rm liq}$ [Pa-s]	$\eta^{\\rm sat}_{\\rm liq}$  [Pa-s]	$\eta^{\\rm sat}_{\\rm liq}$ [Pa-s]\n")
    for rho in range(points):
        sf.write(to_precision(P_t_sat[rho],3)+"	"+to_precision(P_dens_sat[rho],5)+"	")
        if failed_P or rho >= len(P_satmeans):
            sf.write("X	") 
        else:
            sf.write(sig_figs(P_satmeans[rho],P_uncer[rho])+"$_{"+to_precision(P_uncer[rho],2)+"}$	")
        if failed_TA or rho >= len(TA_satmeans):  # We are out of room
            sf.write("X	")
        else:
            sf.write(sig_figs(TA_satmeans[rho],TA_uncer[rho])+"$_{"+to_precision(TA_uncer[rho],2)+"}$	")
        if failed_Tr or rho >= len(Tr_satmeans):  # We are out of room
            sf.write("X\n")
        else:
            sf.write(sig_figs(Tr_satmeans[rho],Tr_uncer[rho])+"$_{"+to_precision(Tr_uncer[rho],2)+"}$\n")


def main():
    
    compound_list = ['C2H6','C3H8','C4H10','C8H18','C12H26','C16H34','C22H46','IC4H10','IC5H12','IC6H14','IC8H18','NEOC5H12','23DMButane','3MPentane']
    compound_nums = {'C2H6':400,'C3H8':400,'C4H10':400,'C8H18':400,'C12H26':400,'C16H34':200,'C22H46':200,'IC4H10':400,'IC5H12':400,'IC6H14':400,'IC8H18':400,'NEOC5H12':400,'23DMButane':400,'3MPentane':400}
    model_list = ['Potoff','TAMie','TraPPE']

    handle_sat_info("C22H46",200)

    for comp in compound_list:
        number = compound_nums[comp]
        handle_sat_info(comp,number)
        handle_T293_info(comp,number)



if __name__ == '__main__':
    '''
    python compare_force_fields.py --comp str --mod str --devPlot --uncer --TDEuncer --Settings
  
    '''

    main()

"""
    comp = "C22H46"
    model = "TraPPE"
    number = 200
    t_sat, dens_sat = load_sim_conditions_sat(comp,model,number)
    t_sim, dens_sim, press_sim, upress_sim = load_sim_conditions_highP(comp,model,number)
    satmeans, satlows, sathighs, satrhos = load_forcefield_data(comp,"Saturation",model,number)
    hpmeans, hplows, hphighs, hprhos = load_forcefield_data(comp,"T293highP",model,number)
    print(t_sat)
    print(dens_sat)
    print(satmeans)
    print(satlows)
    print(sathighs)
    print(satrhos)
    print_T293_info(comp, model, t_sim, dens_sim, press_sim, upress_sim, hpmeans, hplows, hphighs, hprhos)
    print_sat_info(comp, model, t_sat, dens_sat, satmeans, satlows, sathighs, satrhos)
"""


