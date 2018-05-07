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

font = {'size' : '24'}
plt.rc('font',**font)

sim_sat = np.loadtxt('SaturatedSettings.txt',skiprows=1)
Tsat_sim = sim_sat[:,2]

RP_Tsat_eta = np.loadtxt('REFPROP_eta_sat.txt',skiprows=1)
RP_Tsat = RP_Tsat_eta[:,0]
RP_eta_sat = RP_Tsat_eta[:,1]
RP_eta_sim = np.interp(Tsat_sim,RP_Tsat,RP_eta_sat)

TDE_Tsat_eta = np.loadtxt('TDE_eta_sat.txt',skiprows=1)
TDE_Tsat = TDE_Tsat_eta[:,0]
TDE_eta_sat = TDE_Tsat_eta[:,1]

fig,axarr = plt.subplots(ncols=2,nrows=1,figsize=[20,10])

axarr[0].set_yscale("log")
    
for irho in range(5):
    
    try:
        eta_MCMC =  np.loadtxt('GK_eta_boots_rho'+str(irho)) 
        per_error_MCMC = (eta_MCMC-RP_eta_sim[irho])/RP_eta_sim[irho]*100.
        per_error_MCMC_sorted = np.sort(per_error_MCMC)
        i95 = int(0.025*len(per_error_MCMC_sorted))
        per_error_MCMC_plot = per_error_MCMC_sorted[i95:-i95]
        per_error_MCMC_avg = np.mean(per_error_MCMC_plot) #Only take average without the outliers          
        per_error_MCMC_low = per_error_MCMC_avg - per_error_MCMC_plot[0]
        per_error_MCMC_high = per_error_MCMC_plot[-1] - per_error_MCMC_avg
                                
        axarr[0].errorbar(Tsat_sim[irho],np.mean(eta_MCMC),yerr=RP_eta_sim[irho]/100.*np.array([[per_error_MCMC_high,per_error_MCMC_low]]).T,fmt='ro',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
        axarr[1].errorbar(Tsat_sim[irho],per_error_MCMC_avg,yerr=np.array([[per_error_MCMC_high,per_error_MCMC_low]]).T,fmt='ro',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)
    except:
         
        pass

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

fig.savefig('compare_TDE_REFPROP.pdf')

plt.show()