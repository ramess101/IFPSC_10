"""
Compare simulation results with experimental data in TDE and REFPROP correlations
for several force fields and compounds.
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

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--comp",type=str,nargs='+',choices=['C3H8','C4H10','C8H18','IC8H18'],help="Specify the compound(s) to analyze")
    parser.add_argument("-m","--mod",type=str,nargs='+',choices=['Potoff','TraPPE','TAMie','AUA4','Exp6'],help="Specify the model(s) to analyze")
    args = parser.parse_args()
    
    if args.comp and args.mod:
        
        model_list=args.mod
        compound_list=args.comp
        
        #color_list = {'C3H8':'r','C4H10':'b','C8H18':'g','IC8H18':'m'}
        color_list = {'Potoff':'r','TraPPE':'b','TAMie':'g','AUA4':'m'}
        line_list = {'C3H8':'-','C4H10':'--','C8H18':'-.','IC8H18': ':'}
        symbol_list = {'Potoff':'o','TraPPE':'s','TAMie':'^','AUA4':'v'}
        tail_path_list = {'C3H8':'/','C4H10':'_N100/','C8H18':'_N200/','IC8H18':'_N200/'}
        label_list = {'C3H8':r'$n$-C$_3$','C4H10':r'$n$-C$_4$','C8H18':r'$n$-C$_8$','IC8H18':r'$i$-C$_8$'}
        label_position_list = {'C3H8':[330,0.03],'C4H10':[400,0.03],'C8H18':[530,0.03],'IC8H18':[330,0.03]}
        nrhomax = {'C3H8':10,'C4H10':5,'C8H18':5}
        
        eta_avg_all = {0:[],1:[],2:[],3:[],4:[]}
        eta_95_all = {0:[],1:[],2:[],3:[],4:[]}
        
        fig,axarr = plt.subplots(ncols=1,nrows=2,figsize=[10,20])

        axarr[0].set_yscale("log")

        axarr[0].plot([],[],'kx',markersize=10,markeredgewidth=2,mfc='None',label='TDE Exp. Data')
        axarr[0].plot([],[],'k-',linewidth=3,label='REFPROP')

        for model in model_list: axarr[0].plot([],[],color_list[model]+symbol_list[model],markersize=10,markeredgewidth=2,mfc='None',label=model)
        
        for iC, compound in enumerate(compound_list):

            axarr[0].text(label_position_list[compound][0],label_position_list[compound][1],label_list[compound])
            
            for iM, model in enumerate(model_list):
            
                root_path = compound+'/Gromacs/'+model+tail_path_list[compound]
                
                if iM == 0: #Only load the simulation temperatures and the RP/TDE data once
                
                    sim_sat = np.loadtxt(root_path+'SaturatedSettings.txt',skiprows=1)
                    
                    try:
                        Tsat_sim = sim_sat[:,2]
                    except:
                        Tsat_sim = np.array([sim_sat[2]])
                        
                    try:
                        RP_Tsat_eta = np.loadtxt(root_path+'REFPROP_eta_sat.txt',skiprows=1)
                        RP_Tsat = RP_Tsat_eta[:,0]
                        RP_eta_sat = RP_Tsat_eta[:,1]
                        RP_eta_sim = np.interp(Tsat_sim,RP_Tsat,RP_eta_sat)
                        
                        TDE_Tsat_eta = np.loadtxt(root_path+'TDE_eta_sat.txt',skiprows=1)
                        TDE_Tsat = TDE_Tsat_eta[:,0]
                        TDE_eta_sat = TDE_Tsat_eta[:,1]

                        axarr[0].plot(TDE_Tsat,TDE_eta_sat,'kx',markersize=10,markeredgewidth=2,mfc='None')
                        axarr[0].plot(RP_Tsat,RP_eta_sat,'k-',linewidth=3)

                    except:
                        print('Could not find REFPROP or TDE data')
                    
                    nrho = len(Tsat_sim)
                    irho0 = int(nrhomax[compound] - nrho)
    
                for irho in range(nrho):
                    
                    try:
                        eta_MCMC =  np.loadtxt(root_path+'GK_eta_boots_rho'+str(irho+irho0)) 
                        eta_MCMC_sorted = np.sort(eta_MCMC)
                        i95 = int(0.025*len(eta_MCMC_sorted))
                        eta_MCMC_95 = eta_MCMC_sorted[i95:-i95]
                        eta_MCMC_avg = np.mean(eta_MCMC_95) #Only take average without the outliers          
                        eta_MCMC_low = eta_MCMC_avg - eta_MCMC_95[0]
                        eta_MCMC_high = eta_MCMC_95[-1] - eta_MCMC_avg
                                                     
                        #eta_avg_all[irho].append(eta_MCMC_avg) 
                        #eta_95_all[irho].extend(eta_MCMC_95)                                                      

                        per_error_MCMC = (eta_MCMC-RP_eta_sim[irho])/RP_eta_sim[irho]*100.
                        per_error_MCMC_sorted = np.sort(per_error_MCMC)
                        i95 = int(0.025*len(per_error_MCMC_sorted))
                        per_error_MCMC_plot = per_error_MCMC_sorted[i95:-i95]
                        per_error_MCMC_avg = np.mean(per_error_MCMC_plot) #Only take average without the outliers          
                        per_error_MCMC_low = per_error_MCMC_avg - per_error_MCMC_plot[0]
                        per_error_MCMC_high = per_error_MCMC_plot[-1] - per_error_MCMC_avg
                                                
#                        axarr[0].errorbar(Tsat_sim[irho],np.mean(eta_MCMC),yerr=RP_eta_sim[irho]/100.*np.array([[per_error_MCMC_high,per_error_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)    
                        axarr[0].errorbar(Tsat_sim[irho],np.mean(eta_MCMC),yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
                        axarr[1].errorbar(Tsat_sim[irho],per_error_MCMC_avg,yerr=np.array([[per_error_MCMC_high,per_error_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)
                    except:
                        print('Could not find simulation results for '+compound+' '+model+' irho = '+str(irho))
                
        axarr[0].legend()
        
        axarr[0].set_xlabel('Temperature (K)')
        axarr[0].set_ylabel(r'$\eta^{\rm sat}$')
        
        #axarr[0].set_xlim([0.99*Tsat_sim.min(),1.01*Tsat_sim.max()])
        #axarr[0].set_ylim([0.9*RP_eta_sim.min(),1.1*RP_eta_sim.max()])
         
        axarr[1].set_xlabel('Temperature (K)')
        axarr[1].set_ylabel(r'$\left(\eta^{\rm sat}_{\rm sim} - \eta^{\rm sat}_{\rm REFPROP}\right)/\eta^{\rm sat}_{\rm REFPROP} \times 100$%')
        axarr[1].legend()    

        axarr[1].set_xlim(axarr[0].get_xlim())
        
        fig.tight_layout()
        
        fig.savefig('Viscosity_compare_force_fields/compare_force_fields.pdf')
        
        plt.show()  
    
if __name__ == '__main__':
    '''
    python compare_force_fields.py --comp str --mod str
  
    '''

    main()