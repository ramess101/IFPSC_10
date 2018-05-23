"""
Compare simulation results with experimental data in TDE and REFPROP correlations
Must call from "$output_path" with the following files already created:
    SaturatedSettings.txt
    REFPROP_eta_sat.txt
    TDE_eta_sat.txt
    GK_eta_boots_rho"$irho"
"""

from __future__ import division
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse

font = {'size' : '24'}
plt.rc('font',**font)

### Physical constants and conversions
N_A = 6.022140857e23
nm3tom3 = 1e-27
gmtokg = 1e-3

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

Mw_list = {'C8H18':114.23} #[gm/mol]

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--comp",type=str,help="Specify the compound to analyze")
    parser.add_argument("--nrhomax",type=int,help="Specify the maximum number of rho values")
    args = parser.parse_args()

    try:
        comp = args.comp
    except:
        print('Please specify the compound')
        return 0

    Mw = Mw_list[comp]
    
    if args.nrhomax:
        nrhomax = args.nrhomax
    else:
        nrhomax = 5
        
    ### At some point saturated and high pressure should have different names here
    
    sim_sat = np.loadtxt('T347highPSettings.txt',skiprows=1)
    
    try:
        T_sim = sim_sat[:,2]
        L_sim = sim_sat[:,1]
        N_sim = sim_sat[:,0]
    except:
        T_sim = np.array([sim_sat[2]])
        L_sim = np.array([sim_sat[1]])
        N_sim = np.array([sim_sat[0]])

    nrho = nrhomax #len(T_sim)
    
    rho_sim = Lbox_to_rho(L_sim,N_sim,Mw)
    
    RP_T347highP_eta = np.loadtxt('REFPROP_eta_T347highP.txt',skiprows=1)
    RP_press = RP_T347highP_eta[:,1]
    RP_rho = RP_T347highP_eta[:,2]
    RP_eta_T347highP = RP_T347highP_eta[:,-1]
 
    colors_list = {'harmonic':'r','LINCS':'b'}
    shapes_list = {'harmonic':'o','LINCS':'s'}

    fig,axarr = plt.subplots(ncols=2,nrows=1,figsize=[20,10])

    axarr[0].set_yscale("log")
    axarr[1].set_yscale("log")

    axarr[1].plot(RP_press,RP_eta_T347highP,'k-',linewidth=3,label='REFPROP')
    
    ### Read and plot the literature values
    
    eta_TraPPE_Kioupis = np.loadtxt('TraPPE_Kioupis',skiprows=1)
    eta_TraPPE_Nieto = np.loadtxt('TraPPE_Nieto',skiprows=1)

    axarr[0].errorbar(eta_TraPPE_Nieto[:,2],eta_TraPPE_Nieto[:,4],xerr=eta_TraPPE_Nieto[:,3],yerr=eta_TraPPE_Nieto[:,5],fmt='ko',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2,label='Lit. Nieto-Draghi')
    
    axarr[1].plot(eta_TraPPE_Kioupis[:,0],eta_TraPPE_Kioupis[:,1],'ko',mfc='None',markersize=10,markeredgewidth=2,label='Lit. Kioupis')
    axarr[1].errorbar(eta_TraPPE_Nieto[:,0],eta_TraPPE_Nieto[:,4],xerr=eta_TraPPE_Nieto[:,1],yerr=eta_TraPPE_Nieto[:,5],fmt='ks',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)

    axarr[1].errorbar([],[],fmt='ks',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2,label='Lit. Nieto-Draghi')

    irho0 = int(nrhomax - nrho)

    for BondType in ['harmonic', 'LINCS']:

        axarr[1].errorbar([],[],fmt=colors_list[BondType]+shapes_list[BondType],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2,label='This Work, '+BondType)
        
        press_sim = np.loadtxt('TraPPE_N400_'+BondType+'/press_all_log') #[bar]
        press_sim /= 10. #[MPa]
    
        press_sim = press_sim.reshape([nrho,int(len(press_sim)/nrho)])
    
        upress_sim = 1.96*np.std(press_sim,axis=1)
        press_sim = np.mean(press_sim,axis=1)

        for irho in range(nrho):
        
            try:
                eta_MCMC =  np.loadtxt('TraPPE_N400_'+BondType+'/GK_eta_boots_rho'+str(irho+irho0)) 
                eta_MCMC_sorted = np.sort(eta_MCMC)
                i95 = int(0.025*len(eta_MCMC_sorted))
                eta_MCMC_95 = eta_MCMC_sorted[i95:-i95]
                eta_MCMC_avg = np.mean(eta_MCMC_95) #Only take average without the outliers          
                eta_MCMC_low = eta_MCMC_avg - eta_MCMC_95[0]
                eta_MCMC_high = eta_MCMC_95[-1] - eta_MCMC_avg
            
                axarr[0].errorbar(rho_sim[irho],np.mean(eta_MCMC),yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt=colors_list[BondType]+shapes_list[BondType],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
                axarr[1].errorbar(press_sim[irho],np.mean(eta_MCMC),xerr=upress_sim[irho],yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt=colors_list[BondType]+shapes_list[BondType],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
        
            except:
             
                pass
        
    axarr[0].plot([],[],'ro',markersize=10,markeredgewidth=2,mfc='None',label='Simulation')
    axarr[0].plot(RP_rho,RP_eta_T347highP,'k-',linewidth=3,label='REFPROP')
    
    axarr[0].set_xlabel(r'$\rho$ (kg/m$^3$)')
    axarr[0].set_ylabel(r'$\eta$')
     
    axarr[1].legend()

    axarr[1].set_xlabel(r'$P$ (MPa)')
    axarr[1].set_ylabel(r'$\eta$')
     
    axarr[0].set_ylim([0.1,10])
    axarr[1].set_ylim([0.1,10])
    
    fig.tight_layout()

    fig.savefig('compare_TraPPE_lit_T347highP.pdf')

    plt.show()

if __name__ == '__main__':
    '''
    python compare_TraPPE_lit_T347highP.py --comp str --nrhomax int
  
    '''

    main()   