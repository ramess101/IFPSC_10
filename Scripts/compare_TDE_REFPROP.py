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

N_A = 6.022e23
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

Mw_list = {'C3H8':44.096,'C4H10':58.122,'C8H18':114.23}

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--comp",type=str,help="Specify the compound to analyze")
    parser.add_argument("--sat",action='store_true',help="If analyzing saturation conditions")
    parser.add_argument("--T293highP",action='store_true',help="If analyzing T=293 K and high pressures")
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
    
    
    if args.sat and args.T293highP:
        
        print('Must specify either saturation or high pressure conditions')
        return 0
        
    elif args.sat:

        ### At some point saturated and high pressure should have different names here
        
        try:
            sim_sat = np.loadtxt('SaturatedSettings.txt',skiprows=1)
        except:
            sim_sat = np.loadtxt('SaturationSettings.txt',skiprows=1)
        
        try:
            Tsat_sim = sim_sat[:,2]
        except:
            Tsat_sim = np.array([sim_sat[2]])

        nrho = len(Tsat_sim)
        
        RP_Tsat_eta = np.loadtxt('REFPROP_eta_sat.txt',skiprows=1)
        RP_Tsat = RP_Tsat_eta[:,0]
        RP_eta_sat = RP_Tsat_eta[:,1]
        RP_eta_sim = np.interp(Tsat_sim,RP_Tsat,RP_eta_sat)
        
        TDE_Tsat_eta = np.loadtxt('TDE_eta_sat.txt',skiprows=1)
        
        try:
            TDE_Tsat = TDE_Tsat_eta[:,0]
            TDE_eta_sat = TDE_Tsat_eta[:,1]
        except:
            TDE_Tsat = TDE_Tsat_eta[0]
            TDE_eta_sat = TDE_Tsat_eta[1]
        
        fig,axarr = plt.subplots(ncols=2,nrows=1,figsize=[20,10])
    
        axarr[0].set_yscale("log")
        
    elif args.T293highP:
        
        ### At some point saturated and high pressure should have different names here
        
        sim_sat = np.loadtxt('T293highPSettings.txt',skiprows=1)
        
        try:
            T_sim = sim_sat[:,2]
            L_sim = sim_sat[:,1]
            N_sim = sim_sat[:,0]
        except:
            T_sim = np.array([sim_sat[2]])
            L_sim = np.array([sim_sat[1]])
            N_sim = np.array([sim_sat[0]])

        press_sim = np.loadtxt('press_all_log') #[bar]
        press_sim /= 10. #[MPa]

        nrho = len(T_sim)
        
        press_sim = press_sim.reshape([nrho,int(len(press_sim)/nrho)])
        
        upress_sim = 1.96*np.std(press_sim,axis=1)
        press_sim = np.mean(press_sim,axis=1)
        
        rho_sim = Lbox_to_rho(L_sim,N_sim,Mw)
        
        RP_T293highP_eta = np.loadtxt('REFPROP_eta_T293highP.txt',skiprows=1)
        RP_press = RP_T293highP_eta[:,1]
        RP_rho = RP_T293highP_eta[:,2]
        RP_eta_T293highP = RP_T293highP_eta[:,-1]
        RP_eta_rho_sim = np.interp(rho_sim,RP_rho,RP_eta_T293highP)
        RP_eta_press_sim = np.interp(press_sim,RP_press,RP_eta_T293highP)
        
        fig,axarr = plt.subplots(ncols=2,nrows=2,figsize=[20,20])
    
        axarr[0,0].set_yscale("log")
        axarr[0,1].set_yscale("log")
    
    irho0 = int(nrhomax - nrho)
        
    for irho in range(nrho):
        
        try:
            eta_MCMC =  np.loadtxt('GK_eta_boots_rho'+str(irho+irho0)) 
            eta_MCMC_sorted = np.sort(eta_MCMC)
            i95 = int(0.025*len(eta_MCMC_sorted))
            eta_MCMC_95 = eta_MCMC_sorted[i95:-i95]
            eta_MCMC_avg = np.mean(eta_MCMC_95) #Only take average without the outliers          
            eta_MCMC_low = eta_MCMC_avg - eta_MCMC_95[0]
            eta_MCMC_high = eta_MCMC_95[-1] - eta_MCMC_avg
                
            if args.sat:  

                per_error_MCMC = (eta_MCMC-RP_eta_sim[irho])/RP_eta_sim[irho]*100.
                per_error_MCMC_sorted = np.sort(per_error_MCMC)
                i95 = int(0.025*len(per_error_MCMC_sorted))
                per_error_MCMC_plot = per_error_MCMC_sorted[i95:-i95]
                per_error_MCMC_avg = np.mean(per_error_MCMC_plot) #Only take average without the outliers          
                per_error_MCMC_low = per_error_MCMC_avg - per_error_MCMC_plot[0]
                per_error_MCMC_high = per_error_MCMC_plot[-1] - per_error_MCMC_avg                                       
                                                     
                axarr[0].errorbar(Tsat_sim[irho],np.mean(eta_MCMC),yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt='ro',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
                axarr[1].errorbar(Tsat_sim[irho],per_error_MCMC_avg,yerr=np.array([[per_error_MCMC_high,per_error_MCMC_low]]).T,fmt='ro',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)
        
            elif args.T293highP:

                per_error_MCMC = (eta_MCMC-RP_eta_rho_sim[irho])/RP_eta_rho_sim[irho]*100.
                per_error_MCMC_sorted = np.sort(per_error_MCMC)
                i95 = int(0.025*len(per_error_MCMC_sorted))
                per_error_MCMC_plot = per_error_MCMC_sorted[i95:-i95]
                per_error_MCMC_avg = np.mean(per_error_MCMC_plot) #Only take average without the outliers          
                per_error_MCMC_low = per_error_MCMC_avg - per_error_MCMC_plot[0]
                per_error_MCMC_high = per_error_MCMC_plot[-1] - per_error_MCMC_avg 
                
                axarr[0,0].errorbar(rho_sim[irho],np.mean(eta_MCMC),yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt='ro',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
                
                if rho_sim[irho] < RP_rho.max():
                    axarr[1,0].errorbar(rho_sim[irho],per_error_MCMC_avg,yerr=np.array([[per_error_MCMC_high,per_error_MCMC_low]]).T,fmt='ro',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)
                
                per_error_MCMC = (eta_MCMC-RP_eta_press_sim[irho])/RP_eta_press_sim[irho]*100.
                per_error_MCMC_sorted = np.sort(per_error_MCMC)
                i95 = int(0.025*len(per_error_MCMC_sorted))
                per_error_MCMC_plot = per_error_MCMC_sorted[i95:-i95]
                per_error_MCMC_avg = np.mean(per_error_MCMC_plot) #Only take average without the outliers          
                per_error_MCMC_low = per_error_MCMC_avg - per_error_MCMC_plot[0]
                per_error_MCMC_high = per_error_MCMC_plot[-1] - per_error_MCMC_avg 

                axarr[0,1].errorbar(press_sim[irho],np.mean(eta_MCMC),xerr=upress_sim[irho],yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt='ro',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
                
                if press_sim[irho] < RP_press.max():
                    axarr[1,1].errorbar(press_sim[irho],per_error_MCMC_avg,xerr=upress_sim[irho],yerr=np.array([[per_error_MCMC_high,per_error_MCMC_low]]).T,fmt='ro',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)
        
        except:
             
            pass
    
    if args.sat:
    
        axarr[0].plot([],[],'ro',markersize=10,markeredgewidth=2,mfc='None',label='Simulation')
        axarr[0].plot(TDE_Tsat,TDE_eta_sat,'ks',markersize=10,markeredgewidth=2,mfc='None',label='TDE Exp. Data')
        axarr[0].plot(RP_Tsat,RP_eta_sat,'k-',linewidth=3,label='REFPROP')
        axarr[0].legend()
        
        axarr[0].set_xlabel('Temperature (K)')
        axarr[0].set_ylabel(r'$\eta^{\rm sat}$')
        
        #axarr[0].set_xlim([0.99*Tsat_sim.min(),1.01*Tsat_sim.max()])
        #axarr[0].set_ylim([0.9*RP_eta_sim.min(),1.1*RP_eta_sim.max()])
         
        axarr[1].set_xlabel('Temperature (K)')
        axarr[1].set_ylabel(r'$\left(\eta^{\rm sat}_{\rm sim} - \eta^{\rm sat}_{\rm REFPROP}\right)/\eta^{\rm sat}_{\rm REFPROP} \times 100$%')
        axarr[1].legend()   
        
        fig.tight_layout()
    
        fig.savefig('compare_TDE_REFPROP_sat.pdf')
    
        plt.show()
        
    elif args.T293highP:
        
        axarr[0,0].plot([],[],'ro',markersize=10,markeredgewidth=2,mfc='None',label='Simulation')
        axarr[0,0].plot(RP_rho,RP_eta_T293highP,'k-',linewidth=3,label='REFPROP')
        axarr[0,0].legend()
        
        axarr[0,0].set_xlabel(r'$\rho$ (kg/m$^3$)')
        axarr[0,0].set_ylabel(r'$\eta$')
         
        axarr[1,0].set_xlabel(r'$\rho$ (kg/m$^3$)')
        axarr[1,0].set_ylabel(r'$\left(\eta_{\rm sim} - \eta_{\rm REFPROP}\right)/\eta_{\rm REFPROP} \times 100$%')
        #axarr[1,0].legend()   

        axarr[0,1].plot(RP_press,RP_eta_T293highP,'k-',linewidth=3,label='REFPROP')

        axarr[0,1].set_xlabel(r'$P$ (MPa)')
        axarr[0,1].set_ylabel(r'$\eta$')
         
        axarr[1,1].set_xlabel(r'$P$ (MPa)')
        axarr[1,1].set_ylabel(r'$\left(\eta_{\rm sim} - \eta_{\rm REFPROP}\right)/\eta_{\rm REFPROP} \times 100$%')

        axarr[1,0].set_xlim(axarr[0,0].get_xlim())
        axarr[1,1].set_xlim(axarr[0,1].get_xlim())
        
        fig.tight_layout()
    
        fig.savefig('compare_REFPROP_T293highP.pdf')
    
        plt.show()

if __name__ == '__main__':
    '''
    python compare_TDE_REFPROP.py --comp str --nrhomax int --sat --T293highP
  
    '''

    main()   