"""
Plots the dependence of eta on number of molecules
Must call from ~ directory with the following files already created:
    SaturatedSettings.txt
    GK_eta_boots_rho"$irho"
    ~/Viscosity_finite_size_effects/
"""

from __future__ import division
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import linregress
import argparse

font = {'size' : '24'}
plt.rc('font',**font)

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--comp",type=str,help="Specify the compound to analyze")
    parser.add_argument("-m","--mod",type=str,help="Specify the model to analyze")
    parser.add_argument("-N","--Nmol",type=int,nargs='+',help="Specify the number of moles, should be multiple values")
    args = parser.parse_args()
    
    if args.comp and args.mod and args.Nmol:
        
        model=args.mod
        compound=args.comp

        Nmol_list = np.array(args.Nmol)
        directory_list = {400:'/', 100:'_N100/',200:'_N200/'}
        color_list = ['r','b','g','m','c']
        line_list = ['-','--','-','--','-']
        symbol_list = ['o','s','^','v','<','>']
        
        eta_avg_all = {0:[],1:[],2:[],3:[],4:[]}
        eta_95_all = {0:[],1:[],2:[],3:[],4:[]}
        Ninv_all = {0:[],1:[],2:[],3:[],4:[]}
        
        fig,axarr = plt.subplots(ncols=1,nrows=1,figsize=[10,10])
        
        #axarr.set_yscale("log")
        
        for iN in Nmol_list:
            
            root_path = compound+'/Gromacs/'+model+directory_list[iN]
            
            sim_sat = np.loadtxt(root_path+'SaturatedSettings.txt',skiprows=1)
            Tsat_sim = sim_sat[:,2]
            Nmol_sim = sim_sat[0,0]

            nrho = len(Tsat_sim)
            irho0 = int(5 - nrho)
                
            for irho in range(nrho):
                
                try:
                    eta_MCMC =  np.loadtxt(root_path+'GK_eta_boots_rho'+str(irho+irho0)) 
                    eta_MCMC_sorted = np.sort(eta_MCMC)
                    i95 = int(0.025*len(eta_MCMC_sorted))
                    eta_MCMC_95 = eta_MCMC_sorted[i95:-i95]
                    eta_MCMC_avg = np.mean(eta_MCMC_95) #Only take average without the outliers          
                    eta_MCMC_low = eta_MCMC_avg - eta_MCMC_95[0]
                    eta_MCMC_high = eta_MCMC_95[-1] - eta_MCMC_avg
                                                 
                    eta_avg_all[irho].append(eta_MCMC_avg) 
                    eta_95_all[irho].extend(eta_MCMC_95)
                    Ninv_all[irho].extend(Nmol_sim**(-1./3.) for i in range(len(eta_MCMC_95)))                                                         
                    axarr.errorbar(Nmol_sim**(-1./3.),np.mean(eta_MCMC),yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt=color_list[irho]+symbol_list[irho],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
                except:
                     
                    pass
                
                if iN == 400:
        
                    axarr.plot([],[],color_list[irho]+symbol_list[irho]+line_list[irho],linewidth=3,markersize=10,markeredgewidth=2,mfc='None',label=r'$T ='+str(Tsat_sim[irho])+'$ K')
        
        Ninv_plot = np.linspace(0,np.max(Nmol_list**(-1./3.)),1000)
        
        for irho in eta_avg_all:
            
        #    fit = linregress(Nmol_list**(-1./3.),eta_avg_all[irho])
            fit = linregress(Ninv_all[irho],eta_95_all[irho])
            slope = fit[0]
            intercept = fit[1]
            
            plot_fit = intercept + slope*Ninv_plot
            
            axarr.plot(Ninv_plot,plot_fit,color_list[irho]+line_list[irho],linewidth=3)
            axarr.plot([Ninv_plot[0],Ninv_plot[-1]],[np.mean(eta_avg_all[irho]),np.mean(eta_avg_all[irho])],'k:',linewidth=2)
        
        axarr.set_xlabel(r'$N^{-1/3}$')
        axarr.set_ylabel(r'$\eta^{\rm sat}$')
        
        axarr.legend()
        
        #axarr[0].set_xlim([0.99*Tsat_sim.min(),1.01*Tsat_sim.max()])
        #axarr[0].set_ylim([0.9*RP_eta_sim.min(),1.1*RP_eta_sim.max()])
        
        fig.tight_layout()
        
        fig.savefig('Viscosity_finite_size_effects/'+compound+'_'+model+'_finite_size_effects.pdf')
        
        plt.show()

        
    else:
        
        print('Please provide a compound, model, and a list of Nmol')

   
    
if __name__ == '__main__':
    '''
    python plot_finite_size_effects.py --comp str --mod str --Nmol int
  
    '''

    main()   