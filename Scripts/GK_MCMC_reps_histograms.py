# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 07:56:49 2017

@author: ram9
"""

from __future__ import division
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
    
nMCMC = 30    
GK_eta_boots_all = []
GK_low=100.
GK_high=0.

binwidth=0.001

fig, axarr = plt.subplots(1,figsize=(6,6))

for iMCMC in range(nMCMC):
   
   fpath = 'MCMC_'+str(iMCMC)+'/GK_eta_boots'
   
   GK_eta_boots_i = np.loadtxt(fpath)
   GK_eta_boots_i = np.sort(GK_eta_boots_i)
   
   if GK_eta_boots_i[int(97.5*len(GK_eta_boots_i)/100.)] < GK_low:
       
       GK_low = GK_eta_boots_i[int(97.5*len(GK_eta_boots_i)/100.)]
       GK_low_boots = GK_eta_boots_i
       
   if GK_eta_boots_i[int(33.*len(GK_eta_boots_i)/100.)] > GK_high:
       
       GK_high = GK_eta_boots_i[int(33.*len(GK_eta_boots_i)/100.)]
       GK_high_boots = GK_eta_boots_i
   
   GK_eta_boots_all.append(GK_eta_boots_i)                   
   axarr.hist(GK_eta_boots_i,bins=np.arange(min(GK_eta_boots_i),max(GK_eta_boots_i)+binwidth,binwidth),alpha=0.5,normed=True)

GK_eta_boots_all = np.array(GK_eta_boots_all)

GK_eta_boots_all = GK_eta_boots_all.reshape(GK_eta_boots_all.size)
np.savetxt('GK_eta_boots_all',GK_eta_boots_all,fmt='%0.7f')

GK_eta_boots_all = np.loadtxt('GK_eta_boots_all')

GK_eta_boots_all = np.sort(GK_eta_boots_all)

GK_low_95 = GK_eta_boots_all[int(2.5*len(GK_eta_boots_all)/100.)]
GK_high_95 = GK_eta_boots_all[int(97.5*len(GK_eta_boots_all)/100.)]

GK_pu = (GK_high_95-GK_low_95)/(GK_high_95+GK_low_95)*100.

print('Lower 95: '+str(GK_low_95))
print('Upper 95: '+str(GK_high_95))
print('Percent uncertainty in viscosity is: '+str(GK_pu))

hist_max = np.histogram(GK_eta_boots_all,bins=np.arange(min(GK_eta_boots_all),max(GK_eta_boots_all)+binwidth,binwidth))[0].max()
hist_sum = np.histogram(GK_eta_boots_all,bins=np.arange(min(GK_eta_boots_all),max(GK_eta_boots_all)+binwidth,binwidth))[0].sum()

hist_max /= hist_sum
hist_max /= binwidth

GK_eta_MCMC_40_reps = np.loadtxt('../Mie_16_Bayesian_Viscosity_200_MCMC/GK_eta_boots_40reps')
GK_eta_MCMC_200_reps = np.loadtxt('../Mie_16_Bayesian_Viscosity_200_MCMC/GK_eta_boots_200reps')

#axarr[0].hist(GK_eta_boots_all,bins=np.arange(min(GK_eta_boots_all),max(GK_eta_boots_all)+binwidth,binwidth),normed=True)

axarr.hist(GK_high_boots,bins=np.arange(min(GK_high_boots),max(GK_high_boots)+binwidth,binwidth),color='b',normed=True,label='Upper')
axarr.hist(GK_low_boots,bins=np.arange(min(GK_low_boots),max(GK_low_boots)+binwidth,binwidth),color='r',normed=True,label='Lower')

axarr.hist(GK_eta_boots_all,bins=np.arange(min(GK_eta_boots_all),max(GK_eta_boots_all)+binwidth,binwidth),color='k',normed=True,alpha=0.8,label='Combined')
axarr.hist(GK_eta_MCMC_40_reps,bins=np.arange(min(GK_eta_MCMC_40_reps),max(GK_eta_MCMC_40_reps)+binwidth,binwidth),color='g',normed=True,alpha=0.6,label='MCMC 40 subsamples')
axarr.hist(GK_eta_MCMC_200_reps,bins=np.arange(min(GK_eta_MCMC_200_reps),max(GK_eta_MCMC_200_reps)+binwidth,binwidth),color='c',normed=True,alpha=0.6,label='MCMC 200 subsamples')

axarr.plot([GK_low_95,GK_low_95],[0,hist_max],'g:',label='95% credible interval')
axarr.plot([GK_high_95,GK_high_95],[0,hist_max],'g:')

axarr.set_xlabel('Viscosity (cP)')
axarr.set_ylabel('Probability Density (1/cP)')
axarr.legend()
#for ax in axarr:
#
#    ax.set_xlabel('Viscosity (cP)')
#    ax.set_ylabel('Probability Density (1/cP)')
    
plt.tight_layout()

fig.savefig('GK_eta_hist.pdf') 
plt.close()  