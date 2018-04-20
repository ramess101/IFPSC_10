# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Einstein classes

"""

from __future__ import division
import numpy as np 
import os, sys, argparse, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import linregress

TDE_data = np.loadtxt('TDE_viscosity_Ethane.txt',skiprows=1)
T_data = TDE_data[:,0]
eta_data = TDE_data[:,1]

RP_sat = np.loadtxt('RP_eta_sat_Ethane.txt',skiprows=1)
Tsat_RP = RP_sat[:,0]
eta_sat_RP = RP_sat[:,1]

RP_IT = np.loadtxt('RP_eta_IT_Ethane.txt',skiprows=1)
rho_IT_RP = RP_IT[:,0]
eta_IT_RP = RP_IT[:,1]

RP_IT_sim = np.loadtxt('RP_eta_IT_sim_Ethane.txt',skiprows=1)
rho_IT_sim = RP_IT_sim[:,0]
eta_IT_RP_sim = RP_IT_sim[:,1]

TraPPE_data = np.loadtxt('TraPPE_eta_sat.txt',skiprows=1)
Tsat_TraPPE = TraPPE_data[:,0]
eta_sat_TraPPE = TraPPE_data[:,1]
u_low_sat_TraPPE = TraPPE_data[:,2]
u_high_sat_TraPPE = TraPPE_data[:,3]

TraPPE_data = np.loadtxt('TraPPE_eta_IT.txt',skiprows=1)
eta_IT_TraPPE = TraPPE_data[:,0]
u_low_IT_TraPPE = TraPPE_data[:,1]
u_high_IT_TraPPE = TraPPE_data[:,2]

tcut_default = 250

class EinsteinRelationIntegral_MCMC():
    def __init__(self, ilow,ihigh,tcut=tcut_default):
        self.ilow, self.ihigh, self.tcut = ilow, ihigh, tcut
        self.base_path_all = self.gen_base_path()
        self.ERI_MCMC_all = self.gen_ERI_MCMC_all()
        self.nMCMC = len(self.ERI_MCMC_all)
        self.nStates = self.ERI_MCMC_all[0].nStates
        self.nTime = len(self.ERI_MCMC_all[0].ERI_all[0].t_ER)    
        self.ER_MCMC_visc_all, self.EI_MCMC_visc_all, self.t_ERI = self.compile_ERI_MCMC_all()                          
        self.ER_MCMC_visc_avg, self.EI_MCMC_visc_avg, self.EI_MCMC_visc_sig = self.ERI_MCMC_time_avg()
        self.eta_states_ER, self.opt_fit_states = self.calc_eta_states_ER()
        self.eta_states_EI = self.calc_eta_states_EI()
        self.eta_states = self.eta_states_EI
        self.eta_boots, self.eta_low, self.eta_high = self.bootstrap_eta_EI()

    def plot_ERI_all(self):
        """ Plots all the average ER on the same plot """

        self.EI_tcut_MCMC_visc_avg, self.tcut_EI = self.calc_eta_states_EI_t()
        self.EI_t0_MCMC_visc_avg, self.t0_EI = self.calc_eta_states_EI_t(f_t='t0')
        
        nStates, ER_MCMC_visc_avg, t_ER, tcut, EI_tcut_MCMC_visc_avg, tcut_EI, EI_t0_MCMC_visc_avg, t0_EI = self.nStates, self.ER_MCMC_visc_avg, self.t_ERI, self.tcut, self.EI_tcut_MCMC_visc_avg, self.tcut_EI, self.EI_t0_MCMC_visc_avg, self.t0_EI
               
        tplot = np.linspace(0,tcut,1000)
        
        fig = plt.figure(figsize=(6,6))
        
        color_scheme = ['b','r','g','c','m','y','brown','orange','purple','pink','teal','turquoise','violet','indigo','grey','navy','darkgreen','goldenrod','lawngreen','royalblue','tomato', 'maroon']

        for iState in range(nStates):
            
            opt_fit_i = self.opt_fit_states[iState]
            
            ER_data = ER_MCMC_visc_avg[iState,:]
            ER_plot = self.eta_hat(tplot,opt_fit_i)
            
            EI_data_tcut = EI_tcut_MCMC_visc_avg[iState,:]
            EI_data_t0 = EI_t0_MCMC_visc_avg[iState,:]
            
            plt.plot(t_ER,ER_data,'-',color=color_scheme[iState],alpha=0.7)
            plt.plot(tplot,ER_plot,'--',color=color_scheme[iState],alpha=0.7)
            plt.plot(tcut_EI,EI_data_tcut,'-.',color=color_scheme[iState],alpha=0.7)
            plt.plot(t0_EI,EI_data_t0,':',color=color_scheme[iState],alpha=0.7)
            
        plt.xlabel('Time (ps)')
        plt.xlim([0,tcut])
        plt.ylabel('Viscosity (cP)')
        #plt.legend()
                
        fig.savefig('EinsteinRelationIntegral_all.pdf') 
        plt.close()
        
        for iState in range(nStates):
            
            opt_fit_i = self.opt_fit_states[iState]
            
            ER_data = ER_MCMC_visc_avg[iState,:]
            ER_plot = self.eta_hat(tplot,opt_fit_i)
            
            EI_data_tcut = EI_tcut_MCMC_visc_avg[iState,:]
            EI_data_t0 = EI_t0_MCMC_visc_avg[iState,:]

            fig = plt.figure(figsize=(6,6))
            
            plt.plot(t_ER,ER_data,'k-',alpha=0.7,label='Average of slopes')
            plt.plot(tplot,ER_plot,'r--',alpha=0.7,label='Fit to average of slopes')
            plt.plot(tcut_EI,EI_data_tcut,'b-.',alpha=0.7,label='Slope of averages (cut-off)')
            plt.plot(t0_EI,EI_data_t0,'g:',alpha=0.7,label='Slope of averages (t0)')
            
            plt.xlabel('Time (ps)')
            plt.xlim([0,tcut])
            plt.ylabel('Viscosity (cP)')
            plt.legend()
                
            fig.savefig('ERI_state'+str(iState)+'.pdf') 
            plt.close()
        
    def plot_EI_all(self):
        """ Plots all the average ER on the same plot """
        
        nStates, eta_states_EI, tcut, EI_MCMC_visc_avg, t_EI = self.nStates, self.eta_states_EI, self.tcut, self.EI_MCMC_visc_avg, self.t_ERI
        
        tplot = np.linspace(0,tcut,1000)
        
        fig = plt.figure(figsize=(6,6))
        
        color_scheme = ['b','r','g','c','m','y','brown','orange','purple','pink','teal','turquoise','violet','indigo','grey','navy','darkgreen','goldenrod','lawngreen','royalblue','tomato', 'maroon']

        for iState in range(nStates):
            
            eta_state_i = eta_states_EI[iState]
                       
            EI_data = EI_MCMC_visc_avg[iState,:]
            
            EI_plot = eta_state_i*tplot
            
            plt.plot(t_EI,EI_data,'-',color=color_scheme[iState],alpha=0.7)
            plt.plot(tplot,EI_plot,'-.',color=color_scheme[iState],alpha=0.7)
            
        plt.xlabel('Time (ps)')
        plt.xlim([0,tcut])
        plt.ylabel('Einstein Integral (cP)')
                
        fig.savefig('EinsteinIntegral_all.pdf') 
        plt.close()
        
    def plot_sig_all(self):
        """ Plots all the standard deviations on the same plot """
        
        nStates, tcut, EI_MCMC_visc_sig, t_EI = self.nStates, self.tcut, self.EI_MCMC_visc_sig, self.t_ERI
                
        fig = plt.figure(figsize=(6,6))
        
        color_scheme = ['b','r','g','c','m','y','brown','orange','purple','pink','teal','turquoise','violet','indigo','grey','navy','darkgreen','goldenrod','lawngreen','royalblue','tomato', 'maroon']

        for iState in range(nStates):
                                  
            EI_data = EI_MCMC_visc_sig[iState,:]
            
            plt.plot(t_EI,EI_data,'-',color=color_scheme[iState],alpha=0.7)
            
        plt.xlabel('Time (ps)')
        plt.xlim([0,tcut])
        plt.ylabel('Standard Deviation (cP)')
                
        fig.savefig('EinsteinStandardDeviation_all.pdf') 
        plt.close()     
        
    def plot_eta_sat(self):
        
        eta_MCMC, eta_low_MCMC, eta_high_MCMC = self.eta_states, self.eta_low, self.eta_high
        
        Tsat_mask = np.array([False, False, False, False, False, False, False, False, False, True, False, True, False, True, False, True, False, True, False],dtype=bool)
        Tsat = np.array([137.,174.,207.,236.,260.])
        
        eta_sat = eta_MCMC[Tsat_mask]
        eta_low_sat = eta_low_MCMC[Tsat_mask]
        eta_high_sat = eta_high_MCMC[Tsat_mask]
        
        err_low = eta_sat - eta_low_sat
        err_high = eta_high_sat - eta_sat

        fig = plt.figure(figsize=(6,6))
                   
        plt.errorbar(Tsat,eta_sat,yerr=[err_low,err_high],fmt='bs',capsize=2,capthick=2,label='MCMC')
#        plt.plot(Tsat,eta_sat,'bs',label='MCMC')
        plt.errorbar(Tsat_TraPPE,eta_sat_TraPPE,yerr=[u_low_sat_TraPPE,u_high_sat_TraPPE],fmt='r^',capsize=2,capthick=2,label='TraPPE-UA')
        plt.plot(T_data,eta_data,'ko',mfc='None',label='Exp.')
        plt.plot(Tsat_RP,eta_sat_RP,'k--',label='REFPROP')
        plt.ylabel('Viscosity (cP)')
        plt.xlabel('Temperature (K)')
        plt.legend()
        
        fig.savefig('ERI_eta_sat_MCMC.pdf') 
        plt.close()
        
        eta_sat_RP_sim = np.zeros(len(Tsat))
        
        for iSat, Tsat_i in enumerate(Tsat):
            bool_mask = [Tsat_i == Tsat_RP]
            eta_sat_RP_i = eta_sat_RP[bool_mask]
            eta_sat_RP_sim[iSat] = eta_sat_RP_i
                             
        dev_MCMC = (eta_sat-eta_sat_RP_sim)/eta_sat_RP_sim*100.
        dev_TraPPE = (eta_sat_TraPPE-eta_sat_RP_sim)/eta_sat_RP_sim*100.
                     
        pu_low = err_low/eta_sat_RP_sim*100.
        pu_high = err_high/eta_sat_RP_sim*100.   

        pu_low_sat_TraPPE = u_low_sat_TraPPE/eta_sat_RP_sim*100.
        pu_high_sat_TraPPE = u_high_sat_TraPPE/eta_sat_RP_sim*100.          
        
        fig = plt.figure(figsize=(6,6))
                   
        plt.errorbar(Tsat,dev_MCMC,yerr=[pu_low,pu_high],fmt='bs',capsize=2,capthick=2,label='MCMC')
#        plt.plot(Tsat,dev_MCMC,'bs',label='MCMC')
        plt.errorbar(Tsat_TraPPE,dev_TraPPE,yerr=[pu_low_sat_TraPPE,pu_high_sat_TraPPE],fmt='r^',capsize=2,capthick=2,label='TraPPE-UA')
        #plt.plot(T_data,eta_data,'ko',mfc='None',label='Exp.')
        #plt.plot(Tsat_RP,eta_sat_RP,'k--',label='REFPROP')
        plt.ylabel('Percent Deviation in Viscosity')
        plt.xlabel('Temperature (K)')
        plt.legend()
        
        fig.savefig('ERI_dev_eta_sat_MCMC.pdf') 
        plt.close()
        
    def plot_eta_IT(self):
        
        eta_MCMC, eta_low_MCMC, eta_high_MCMC = self.eta_states, self.eta_low, self.eta_high
        
        IT_mask = np.array([True, True, True, True, True, True, True, True, True, False, False, False, False, False, False, False, False, False, False],dtype=bool)
        
        eta_IT = eta_MCMC[IT_mask]
        eta_low_IT = eta_low_MCMC[IT_mask]
        eta_high_IT = eta_high_MCMC[IT_mask]
#        
        err_low = eta_IT - eta_low_IT
        err_high = eta_high_IT - eta_IT

        fig = plt.figure(figsize=(6,6))
                   
        plt.errorbar(rho_IT_sim,eta_IT,yerr=[err_low,err_high],fmt='bs',capsize=2,capthick=2,label='MCMC')
#        plt.plot(rho_IT_sim,eta_IT,'bs',label='MCMC')
        plt.errorbar(rho_IT_sim,eta_IT_TraPPE,yerr=[u_low_IT_TraPPE,u_high_IT_TraPPE],fmt='r^',capsize=2,capthick=2,label='TraPPE-UA')
        plt.plot(rho_IT_RP,eta_IT_RP,'k--',label='REFPROP')
        plt.ylabel('Viscosity (cP)')
        plt.xlabel('Temperature (K)')
        plt.legend()
        
        fig.savefig('ERI_eta_IT_MCMC.pdf') 
        plt.close()
                                     
        dev_MCMC = (eta_IT-eta_IT_RP_sim)/eta_IT_RP_sim*100.
        dev_TraPPE = (eta_IT_TraPPE-eta_IT_RP_sim)/eta_IT_RP_sim*100.
                     
        pu_low = err_low/eta_IT_RP_sim*100.
        pu_high = err_high/eta_IT_RP_sim*100.   

        pu_low_IT_TraPPE = u_low_IT_TraPPE/eta_IT_RP_sim*100.
        pu_high_IT_TraPPE = u_high_IT_TraPPE/eta_IT_RP_sim*100.          
        
        fig = plt.figure(figsize=(6,6))
                   
        plt.errorbar(rho_IT_sim,dev_MCMC,yerr=[pu_low,pu_high],fmt='bs',capsize=2,capthick=2,label='MCMC')
#        plt.plot(rho_IT_sim,dev_MCMC,'bs',label='MCMC')
        plt.errorbar(rho_IT_sim,dev_TraPPE,yerr=[pu_low_IT_TraPPE,pu_high_IT_TraPPE],fmt='r^',capsize=2,capthick=2,label='TraPPE-UA')
        plt.ylabel('Percent Deviation in Viscosity')
        plt.xlabel('Temperature (K)')
        plt.legend()
        
        fig.savefig('ERI_dev_eta_IT_MCMC.pdf') 
        plt.close()
            

    def compile_ERI_MCMC_all(self):
        """ Compile the Green-Kubo values for all time, MCMC samples, and states """
        nStates,nTime,nMCMC = self.nStates, self.nTime, self.nMCMC
        
        ERI_MCMC_all = self.ERI_MCMC_all
        
        ER_MCMC_visc_all = np.zeros([nStates,nTime,nMCMC])
        EI_MCMC_visc_all = np.zeros([nStates,nTime,nMCMC])
        
        for iState in range(nStates):
            
            for iMCMC in range(nMCMC):
                
                ER_MCMC_visc_all[iState,:,iMCMC] = ERI_MCMC_all[iMCMC].ERI_all[iState].visc_ER
                EI_MCMC_visc_all[iState,:,iMCMC] = ERI_MCMC_all[iMCMC].ERI_all[iState].visc_EI
                                
        t_ERI = ERI_MCMC_all[0].ERI_all[0].t_ER
            
        return ER_MCMC_visc_all, EI_MCMC_visc_all, t_ERI
    
    def ERI_MCMC_time_avg(self):
        """ Averages the viscosity at each time for different MCMC samples """
        ER_MCMC_visc_all, EI_MCMC_visc_all = self.ER_MCMC_visc_all, self.EI_MCMC_visc_all
        
        ER_MCMC_visc_avg = np.mean(ER_MCMC_visc_all,axis=2)
        EI_MCMC_visc_avg = np.mean(EI_MCMC_visc_all,axis=2)
        EI_MCMC_visc_sig = np.std(EI_MCMC_visc_all,axis=2)

        np.savetxt('ER_MCMC_avg',ER_MCMC_visc_avg.T,fmt='%0.7f')
        np.savetxt('EI_MCMC_avg',EI_MCMC_visc_avg.T,fmt='%0.7f')
        np.savetxt('EI_MCMC_sig',EI_MCMC_visc_avg.T,fmt='%0.7f')
        
        return ER_MCMC_visc_avg, EI_MCMC_visc_avg,EI_MCMC_visc_sig            
    
    def calc_eta_states_ER(self):
        
        ER_MCMC_visc_avg,nStates,t_ER = self.ER_MCMC_visc_avg,self.nStates,self.t_ERI
      
        eta_states = np.zeros(nStates)
        opt_fit_states = np.zeros([nStates,4])
        
        for iState in range(nStates):
            
            eta_avg = ER_MCMC_visc_avg[iState,:]
            opt_fit = self.fit_eta(t_ER,eta_avg)
            eta_states[iState] = self.calc_eta_estimate(opt_fit)
            opt_fit_states[iState,:] = opt_fit
        
        np.savetxt('ER_eta_states',eta_states,fmt='%0.7f')
        
        return eta_states, opt_fit_states
    
    def calc_eta_states_EI_t(self,f_t='tcut'):
        
        EI_MCMC_visc_avg,nStates,t_EI = self.EI_MCMC_visc_avg,self.nStates,self.t_ERI
              
        n_t = len(t_EI)
        
        if f_t == 'tcut':
            
            t_cuts = np.linspace(1,250,50)
            
        elif f_t == 't0':
             
            t_cuts = np.linspace(0,249,50)
        
        n_t = len(t_cuts)
        
        EI_t = np.zeros([nStates,n_t])
        
        for iState in range(nStates):
            
            EI_avg = EI_MCMC_visc_avg[iState,:]
            
            for i_t, tcut in enumerate(t_cuts):
            
                if f_t == 'tcut':
                
                    EI_data = EI_avg[t_EI<tcut]
                    t_data = t_EI[t_EI<tcut]
                
                elif f_t == 't0':
                
                    EI_data = EI_avg[t_EI>tcut]
                    t_data = t_EI[t_EI>tcut]
                
#                slope_fit = np.polyfit(t_data,EI_data,1)[0]
                slope_fit = linregress(t_data,EI_data)[0]
                EI_t[iState,i_t] = slope_fit
        
        return EI_t, t_cuts
    
    def calc_eta_states_EI(self):
        
        EI_MCMC_visc_avg,nStates,t_EI = self.EI_MCMC_visc_avg,self.nStates,self.t_ERI
      
        eta_states = np.zeros(nStates)
        
        for iState in range(nStates):
            
            EI_avg = EI_MCMC_visc_avg[iState,:]
#            slope_fit = np.polyfit(t_EI,EI_avg,1)
            slope_fit = linregress(t_EI,EI_avg)[0]
            eta_states[iState] = slope_fit
        
        np.savetxt('EI_eta_states',eta_states,fmt='%0.7f')
        
        return eta_states
    
    def bootstrap_eta_EI(self):
        
        EI_MCMC_visc_avg,nStates,t_EI = self.EI_MCMC_visc_avg,self.nStates,self.t_ERI
        
        nBoots = 200
      
        eta_boots = np.zeros([nStates,nBoots])
        eta_low = np.zeros(nStates)
        eta_high = np.zeros(nStates)
        
        tcut_range = np.linspace(50.,250.,10000)
        t0_range = np.linspace(0,200,10000)
        
        for iState in range(nStates):
            
            EI_avg = EI_MCMC_visc_avg[iState,:]
            
            for iBoots in range(nBoots):
                
                if np.random.random()<0.5:
                
                    tcut = np.random.choice(tcut_range)
                    
                    EI_data = EI_avg[t_EI<tcut]
                    t_data = t_EI[t_EI<tcut]
                
                else:
                    
                    tcut = np.random.choice(t0_range)
                    
                    EI_data = EI_avg[t_EI>tcut]
                    t_data = t_EI[t_EI>tcut]
                
                slope_fit = linregress(t_data,EI_data)[0]
                eta_boots[iState,iBoots] = slope_fit       
            
            eta_boots_sorted = np.sort(np.array(eta_boots[iState,:]))
            
            ilow = int(nBoots*2.5/100.)
            ihigh = int(nBoots*97.5/100.)
            
            eta_low[iState] = eta_boots_sorted[ilow]
            eta_high[iState] = eta_boots_sorted[ihigh]
            
        np.savetxt('EI_eta_boots',eta_boots[:,:200],fmt='%0.7f')
            
        return eta_boots, eta_low, eta_high
    
    def plot_bootstraps(self):
        """ Plots the bootstrapped histogram """
        eta_boots_MCMC, eta_low_MCMC, eta_high_MCMC, eta_states = self.eta_boots, self.eta_low, self.eta_high, self.eta_states
        
        nbins = int(len(eta_boots_MCMC.T)/5)
        
        for iState in range(self.nStates):
            
            eta_boots = eta_boots_MCMC[iState,:]
            eta_low = eta_low_MCMC[iState]
            eta_high = eta_high_MCMC[iState]
            eta_state = eta_states[iState]
        
            hist_bins = np.histogram(eta_boots,bins=nbins)[1]
            dbin = hist_bins[1]-hist_bins[0]
            hist_max = np.histogram(eta_boots,bins=nbins)[0].max()
            hist_sum = np.histogram(eta_boots,bins=nbins)[0].sum()
    
            hist_max /= hist_sum
            hist_max /= dbin
    
            fig = plt.figure(figsize=(6,6))
                   
            plt.hist(eta_boots,bins=nbins,normed=True,color='b')
            plt.plot([eta_low,eta_low],[0,hist_max],'k--',label=r'95% confidence interval')
            plt.plot([eta_high,eta_high],[0,hist_max],'k--')
            plt.plot([eta_state,eta_state],[0,hist_max],'r-',label='Estimated Viscosity')
            plt.xlabel('Viscosity (cP)')
            plt.ylabel('Probability Density (1/cP)')
            plt.legend()
            
            fig.savefig('ER_bootstraps_state'+str(iState)+'.pdf') 
            plt.close()
    
    def eta_hat(self,t,params):
        """ Model suggested by Maginn for fitting viscosity """
        
        A = params[0]
        alpha = params[1]
        tau1 = params[2]
        tau2 = params[3]
        
        eta_t = A*alpha*tau1*(1-np.exp(-t/tau1))+A*(1-alpha)*tau2*(1-np.exp(-t/tau2))
        return eta_t
    
    def fit_eta(self,t_data,eta_data,tcut=tcut_default,w8_data=None):
        """ Fits the viscosity data to correlation with assigned weights """
        
        eta_data = eta_data[t_data<tcut]
        t_data = t_data[t_data<tcut]
        
        if w8_data == None:
            w8_data = t_data**(-0.5)
        
        eta_t = lambda params: self.eta_hat(t_data,params)
        
        SSE = lambda params: np.sum(((eta_t(params) - eta_data)*w8_data)**2)
                             
        guess = [2, 0.1, 0.5, 0.1]
        
        bnds=((0,None),(0,None),(1e-5,None),(1e-5,None))
                     
        opt = minimize(SSE,guess,bounds=bnds)

        opt_fit = opt.x                     
    
        return opt_fit
    
    def calc_eta_estimate(self,opt_fit):
        
        A = opt_fit[0]
        alpha = opt_fit[1]
        tau1 = opt_fit[2]
        tau2 = opt_fit[3]
        
        eta_estimate = A*alpha*tau1+A*(1-alpha)*tau2
        return eta_estimate   
    
    def gen_base_path(self):
        ilow, ihigh = self.ilow, self.ihigh
        
        base_path_all = []
        
        for iMCMC in np.arange(ilow,ihigh+1):
            
            base_path_all.append('MCMC_'+str(iMCMC)+'/')
   
        return base_path_all
           
    def gen_ERI_MCMC_all(self):
        
        base_path_all = self.base_path_all
        ERI_MCMC_all = []
        
        for base_path_i in base_path_all:
            
            ERI_MCMC_all.append(EinsteinRelationIntegral_states(base_path=base_path_i))
            
        return ERI_MCMC_all
       

class EinsteinRelationIntegral_states():
    def __init__(self,base_path='',fpath_all=None):
        try:
            self.fpath_all = self.gen_fpath_all()
        except:
            self.fpath_all = fpath_all
        self.base_path = base_path
        self.total_path_all = [self.base_path + fpath for fpath in self.fpath_all]    
        self.nStates = len(self.fpath_all)
        self.ERI_all = self.gen_ERI_all()
        
    def gen_fpath_all(self):
        #Read in the simulation specifications
    
        ITIC = np.array(['Isotherm', 'Isochore'])
          
        Temp_ITIC = {'Isochore':[],'Isotherm':[]}
        rho_ITIC = {'Isochore':[],'Isotherm':[]}
        Nmol = {'Isochore':[],'Isotherm':[]}
        Temps = {'Isochore':[],'Isotherm':[]}
        rhos_ITIC = {'Isochore':[],'Isotherm':[]}
        nTemps = {'Isochore':[],'Isotherm':[]}
        nrhos = {'Isochore':[],'Isotherm':[]}
        
        for run_type in ITIC:
        
            run_type_Settings = np.loadtxt(run_type+'Settings.txt',skiprows=1)
            
            Nmol[run_type] = run_type_Settings[:,0]
            Lbox = run_type_Settings[:,1] #[nm]
            Temp_ITIC[run_type] = run_type_Settings[:,2] #[K]
            Vol = Lbox**3 #[nm3]
            rho_ITIC[run_type] = Nmol[run_type] / Vol #[molecules/nm3]
            rhos_ITIC[run_type] = np.unique(rho_ITIC[run_type])
        
            nrhos[run_type] = len(rhos_ITIC[run_type])
            Temps[run_type] = np.unique(Temp_ITIC[run_type])
            nTemps[run_type] = len(Temps[run_type]) 
             
        nTemps['Isochore']=2 #Need to figure out how to get this without hardcoding
        
        # Create a list of all the file paths (without the reference directory, just the run_type, rho, Temp)
        
        fpath_all = []
        
        for run_type in ITIC: 
        
            for irho  in np.arange(0,nrhos[run_type]):
        
                for iTemp in np.arange(0,nTemps[run_type]):
        
                    if run_type == 'Isochore':
        
                        fpath_all.append(run_type+'/rho'+str(irho)+'/T'+str(iTemp)+'/NVT_eq/NVT_prod/NVT_vis/')
        
                    else:
        
                        fpath_all.append(run_type+'/rho_'+str(irho)+'/NVT_eq/NVT_prod/NVT_vis/')
                        
        return fpath_all
                
    def gen_ERI_all(self):
        
        total_path_all = self.total_path_all
        ERI_all = []
        
        for path in total_path_all:
            
            ERI_all.append(EinsteinRelationIntegral(path))
    
        return ERI_all
    

class EinsteinRelationIntegral():
    def __init__(self, fpath):
        self.fpath = fpath
        self.t_ER, self.visc_ER = self.load_visc_ER()
        self.visc_EI = self.load_visc_EI()
        
    def load_visc_ER(self):
        """ Reads in the viscosities for the Einstein relation """
        fpath = self.fpath
        g_start = 17
        t_visc_ER = open(fpath+'evisco.xvg','r').readlines()[g_start:] #Read all lines starting at g_start
        n_t = len(t_visc_ER)
                    
        t_ER = np.zeros(n_t)
        visc_ER = np.zeros(n_t)
         
        for i_t in range(n_t):
            t_ER[i_t] = float(t_visc_ER[i_t].split()[0])
            visc_ER[i_t] = float(t_visc_ER[i_t].split()[4])
            
        visc_ER *= 1000. #Convert from [kg/m/s] to [cP]
        return t_ER, visc_ER
    
    def load_visc_EI(self):
        """ Reads in the viscosities for the Einstein integral """
        fpath = self.fpath
        g_start = 17
        t_visc_EI = open(fpath+'eviscoi.xvg','r').readlines()[g_start:] #Read all lines starting at g_start
        n_t = len(t_visc_EI)
                    
        t_EI = np.zeros(n_t)
        visc_EI = np.zeros(n_t)
         
        for i_t in range(n_t):
            t_EI[i_t] = float(t_visc_EI[i_t].split()[0])
            visc_EI[i_t] = float(t_visc_EI[i_t].split()[4])
            
        visc_EI *= 1000. #Convert from [kg/m/s] to [cP]
        return visc_EI

#fpath = 'TraPPE_IT_rho8/'   
#TraPPE_IT_rho8 = EinsteinRelation(fpath)

#plt.plot(TraPPE_IT_rho8.t_ER,TraPPE_IT_rho8.visc_ER)
#plt.show()
#print(TraPPE_IT_rho8.t_ER[0])

#bpath = 'MCMC_0/'
#MCMC_0 = EinsteinRelation_states(base_path=bpath)
    
#ERI_MCMC = EinsteinRelationIntegral_MCMC(0,1)
#print(ERI_MCMC.eta_states_EI)
#print(ERI_MCMC.eta_states_ER)
#ERI_MCMC.plot_ERI_all()
#ER_MCMC.plot_eta_IT()
#ER_MCMC.plot_bootstraps()
#print(ER_MCMC.eta_states)
#for iState in range(ER_MCMC.nStates):
#    
#    plt.plot(ER_MCMC.t_ER,ER_MCMC.ER_MCMC_visc_avg[iState,:])
#    
#plt.show()
    
def main():
    
#    fpath = 'TraPPE_IT_rho8/'
#    TraPPE_IT_rho8 = EinsteinRelation(fpath)
#    
#    plt.plot(TraPPE_IT_rho8.t_ER,TraPPE_IT_rho8.visc_ER)
#    plt.show()
#
#fpath = ['Potoff_IT_rho8/','TraPPE_IT_rho8/']
#multi_ER = EinsteinRelation_states(fpath_all=fpath)
#
#plt.plot(multi_ER.ER_all[0].t_ER,multi_ER.ER_all[0].visc_ER)
#plt.plot(multi_ER.ER_all[1].t_ER,multi_ER.ER_all[1].visc_ER)
#plt.show()

#    ER_MCMC = EinsteinRelation_MCMC(0,1)
#    ER_MCMC.plot_eta_IT()
#    ER_MCMC.plot_bootstraps()
#    print(ER_MCMC.eta_states)
#    for iState in range(ER_MCMC.nStates):
#        
#        plt.plot(ER_MCMC.t_ER,ER_MCMC.ER_MCMC_visc_avg[iState,:])
#        
#    plt.show()
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-il","--ilow",type=int,help="Specify the lowest iMCMC value")
    parser.add_argument("-ih","--ihigh",type=int,help="Specify the high iMCMC value")
    args = parser.parse_args()
    
    MCMC_Mie = EinsteinRelationIntegral_MCMC(args.ilow,args.ihigh)
    MCMC_Mie.plot_ERI_all()
    MCMC_Mie.plot_EI_all()
    MCMC_Mie.plot_sig_all()
    MCMC_Mie.plot_eta_sat()
    MCMC_Mie.plot_eta_IT()
    MCMC_Mie.plot_bootstraps()
#    
if __name__ == '__main__':
    '''
    python Einstein_analyze.py --ilow int --ihigh int
  
    '''

    main()   
