# -*- coding: utf-8 -*-
"""
Green-Kubo classes

"""

from __future__ import division
import numpy as np 
import os, sys, argparse, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize

try:
    
    TDE_data = np.loadtxt('TDE_eta_sat.txt',skiprows=1)
    T_data = TDE_data[:,0]
    eta_data = TDE_data[:,1]
    
    RP_sat = np.loadtxt('REFPROP_eta_sat.txt',skiprows=1)
    Tsat_RP = RP_sat[:,0]
    eta_sat_RP = RP_sat[:,1]
    
#    RP_IT = np.loadtxt('REFPROP_eta_IT.txt',skiprows=1)
#    rho_IT_RP = RP_IT[:,0]
#    eta_IT_RP = RP_IT[:,1]
    
#    RP_IT_sim = np.loadtxt('REFPROP_eta_IT_sim.txt',skiprows=1)
#    rho_IT_sim = RP_IT_sim[:,0]
#    eta_IT_RP_sim = RP_IT_sim[:,1]
    
except:
    
    pass

sim_sat = np.loadtxt('SaturatedSettings.txt',skiprows=1)
Tsat = sim_sat[:,2]

tcut_default = 10
tlow_default = 2

class GreenKubo_MCMC():
    def __init__(self, ilow,ihigh,tcut=tcut_default):
        self.ilow, self.ihigh, self.tcut = ilow, ihigh, tcut
        self.base_path_all = self.gen_base_path()
        self.GK_MCMC_all = self.gen_GK_MCMC_all()
        self.nMCMC = len(self.GK_MCMC_all)
        self.nStates = self.GK_MCMC_all[0].nStates
        self.nTime = len(self.GK_MCMC_all[0].GK_all[0].t_GK)    
        self.GK_MCMC_visc_all, self.t_GK = self.compile_GK_MCMC_all()     
        self.w8_model = self.w8_hat()                     
        self.GK_MCMC_visc_avg = self.GK_MCMC_time_avg()
        self.eta_states, self.opt_fit_states = self.calc_eta_states()
        self.eta_boots, self.eta_low, self.eta_high = self.bootstrap_eta()
        
    def w8_hat(self,fit=False):
        t_GK = self.t_GK
        if not fit:
            w8_model = t_GK**(-0.5)
        return w8_model

    def plot_GK_all(self):
        """ Plots all the average GK on the same plot """
        
        nStates, GK_MCMC_visc_avg, t_GK, tcut = self.nStates, self.GK_MCMC_visc_avg, self.t_GK, self.tcut
        
        tplot = np.linspace(0,tcut,1000)
        
        fig = plt.figure(figsize=(6,6))
        
        color_scheme = ['b','r','g','c','m','y','brown','orange','purple','pink','teal','turquoise','violet','indigo','grey','navy','darkgreen','goldenrod','lawngreen','royalblue','tomato', 'maroon']

        for iState in range(nStates):
            
            opt_fit_i = self.opt_fit_states[iState]
            
            GK_data = GK_MCMC_visc_avg[iState,:]
            GK_plot = self.eta_hat(tplot,opt_fit_i)
            
 #           plt.plot(t_GK[t_GK<tcut],GK_data[t_GK<tcut],'-',color=color_scheme[iState],alpha=0.7)
            plt.plot(t_GK,GK_data,'-',color=color_scheme[iState],alpha=0.7)
            plt.plot(tplot,GK_plot,'--',color=color_scheme[iState+1],alpha=0.7)
            
        plt.xlabel('Time (ps)')
        plt.xlim([0,tcut])
        plt.ylabel('Viscosity (cP)')
        #plt.legend()
                
        fig.savefig('GK_all.pdf') 
        plt.close() 
        
    def plot_eta_sat(self):
        
        eta_MCMC, eta_low_MCMC, eta_high_MCMC = self.eta_states, self.eta_low, self.eta_high
        
        Tsat_mask = np.array([False, False, False, False, False, False, False, False, False, True, False, True, False, True, False, True, False, True, False],dtype=bool)
        
        eta_sat = eta_MCMC[Tsat_mask]
        eta_low_sat = eta_low_MCMC[Tsat_mask]
        eta_high_sat = eta_high_MCMC[Tsat_mask]
        
        err_low = eta_sat - eta_low_sat
        err_high = eta_high_sat - eta_sat

        fig = plt.figure(figsize=(6,6))
                   
        plt.errorbar(Tsat,eta_sat,yerr=[err_low,err_high],fmt='bs',capsize=2,capthick=2,label='MCMC')
#        plt.plot(Tsat,eta_sat,'bs',label='MCMC')
        plt.plot(T_data,eta_data,'ko',mfc='None',label='Exp.')
        plt.plot(Tsat_RP,eta_sat_RP,'k--',label='REFPROP')
        plt.ylabel('Viscosity (cP)')
        plt.xlabel('Temperature (K)')
        plt.legend()
        
        fig.savefig('GK_eta_sat_MCMC.pdf') 
        plt.close()
        
        eta_sat_RP_sim = np.zeros(len(Tsat))
        
        for iSat, Tsat_i in enumerate(Tsat):
            bool_mask = [Tsat_i == Tsat_RP]
            eta_sat_RP_i = eta_sat_RP[bool_mask]
            eta_sat_RP_sim[iSat] = eta_sat_RP_i
                             
        dev_MCMC = (eta_sat-eta_sat_RP_sim)/eta_sat_RP_sim*100.
                     
        pu_low = err_low/eta_sat_RP_sim*100.
        pu_high = err_high/eta_sat_RP_sim*100.           
        
        fig = plt.figure(figsize=(6,6))
                   
        plt.errorbar(Tsat,dev_MCMC,yerr=[pu_low,pu_high],fmt='bs',capsize=2,capthick=2,label='MCMC')
#        plt.plot(Tsat,dev_MCMC,'bs',label='MCMC')
        #plt.plot(T_data,eta_data,'ko',mfc='None',label='Exp.')
        #plt.plot(Tsat_RP,eta_sat_RP,'k--',label='REFPROP')
        plt.ylabel('Percent Deviation in Viscosity')
        plt.xlabel('Temperature (K)')
        plt.legend()
        
        fig.savefig('GK_dev_eta_sat_MCMC.pdf') 
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
        plt.plot(rho_IT_RP,eta_IT_RP,'k--',label='REFPROP')
        plt.ylabel('Viscosity (cP)')
        plt.xlabel('Temperature (K)')
        plt.legend()
        
        fig.savefig('GK_eta_IT_MCMC.pdf') 
        plt.close()
                                     
        dev_MCMC = (eta_IT-eta_IT_RP_sim)/eta_IT_RP_sim*100.
                     
        pu_low = err_low/eta_IT_RP_sim*100.
        pu_high = err_high/eta_IT_RP_sim*100.         
        
        fig = plt.figure(figsize=(6,6))
                   
        plt.errorbar(rho_IT_sim,dev_MCMC,yerr=[pu_low,pu_high],fmt='bs',capsize=2,capthick=2,label='MCMC')
#        plt.plot(rho_IT_sim,dev_MCMC,'bs',label='MCMC')
        plt.ylabel('Percent Deviation in Viscosity')
        plt.xlabel('Temperature (K)')
        plt.legend()
        
        fig.savefig('GK_dev_eta_IT_MCMC.pdf') 
        plt.close()
            
    def compile_GK_MCMC_all(self):
        """ Compile the Green-Kubo values for all time, MCMC samples, and states """
        nStates,nTime,nMCMC = self.nStates, self.nTime, self.nMCMC
        
        GK_MCMC_all = self.GK_MCMC_all
        
        GK_MCMC_visc_all = np.zeros([nStates,nTime,nMCMC])
        
        for iState in range(nStates):
            
            for iMCMC in range(nMCMC):
                
                GK_MCMC_visc_all[iState,:,iMCMC] = GK_MCMC_all[iMCMC].GK_all[iState].visc_GK
                                
        t_GK = GK_MCMC_all[0].GK_all[0].t_GK
            
        return GK_MCMC_visc_all, t_GK
    
    def GK_MCMC_time_avg(self):
        """ Averages the viscosity at each time for different MCMC samples """
        GK_MCMC_visc_all = self.GK_MCMC_visc_all
        
        GK_MCMC_visc_avg = np.mean(GK_MCMC_visc_all,axis=2)

        np.savetxt('GK_MCMC_avg',GK_MCMC_visc_avg.T,fmt='%0.7f')
        
        return GK_MCMC_visc_avg            
    
    def calc_eta_states(self):
        
        GK_MCMC_visc_avg,nStates,t_GK, w8_model = self.GK_MCMC_visc_avg,self.nStates,self.t_GK,self.w8_model
      
        eta_states = np.zeros(nStates)
        opt_fit_states = np.zeros([nStates,4])
                
        for iState in range(nStates):
            
            eta_avg = GK_MCMC_visc_avg[iState,:]
            opt_fit = self.fit_eta(t_GK,eta_avg,w8_model)
            eta_states[iState] = self.calc_eta_estimate(opt_fit)
            opt_fit_states[iState,:] = opt_fit
        
        np.savetxt('GK_eta_states',eta_states,fmt='%0.7f')
        
        return eta_states, opt_fit_states
    
    def bootstrap_eta(self):
        
        GK_MCMC_visc_avg,nStates,t_GK, w8_model,GK_MCMC_visc_all = self.GK_MCMC_visc_avg,self.nStates,self.t_GK,self.w8_model, self.GK_MCMC_visc_all
        
        nBoots = 2
      
        eta_boots = np.zeros([nStates,nBoots])
        eta_low = np.zeros(nStates)
        eta_high = np.zeros(nStates)
        
        tcut_range = np.linspace(0.8*self.tcut,10*self.tcut)
        b_low = 0.4
        b_high = 0.6
        
        for iState in range(nStates):
            
            eta_avg = GK_MCMC_visc_avg[iState,:]
            eta_state_all = GK_MCMC_visc_all[iState,:,:]
            
            for iBoots in range(nBoots):
                
#                w8_boots = np.random.random(self.nMCMC) #This has a random weight for the different runs, essentially it randomly samples which runs to include 
#                w8_boots = np.random.randint(0,2,self.nMCMC) # The problem with this method is that you might have all 0's
                w8_boots = np.random.randint(0,self.nMCMC,self.nMCMC) # This method tries to ensure that a large number of runs are included
                while np.sum(w8_boots) == 0: w8_boots = np.random.randint(0,2,self.nMCMC)
 
                eta_avg = np.sum(eta_state_all*w8_boots,axis=1)/np.sum(w8_boots)
                
                tcut = np.random.choice(tcut_range)
                b_random = np.random.uniform(b_low,b_high)
                w8_model = t_GK**(-b_random)
            
                opt_fit = self.fit_eta(t_GK,eta_avg,w8_model,tcut)
                eta_boots[iState,iBoots] = self.calc_eta_estimate(opt_fit)
            
            eta_boots_sorted = np.sort(np.array(eta_boots[iState,:]))
            
            ilow = int(nBoots*2.5/100.)
            ihigh = int(nBoots*97.5/100.)
            
            eta_low[iState] = eta_boots_sorted[ilow]
            eta_high[iState] = eta_boots_sorted[ihigh]
            
        np.savetxt('GK_eta_boots',eta_boots[:,:200],fmt='%0.7f')
            
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
            
            fig.savefig('GK_bootstraps_state'+str(iState)+'.pdf') 
            plt.close()
    
    def eta_hat(self,t,params):
        """ Model suggested by Maginn for fitting viscosity """
        
        A = params[0]
        alpha = params[1]
        tau1 = params[2]
        tau2 = params[3]
        
        eta_t = A*alpha*tau1*(1-np.exp(-t/tau1))+A*(1-alpha)*tau2*(1-np.exp(-t/tau2))
        return eta_t
    
    def fit_eta(self,t_data,eta_data,w8_data,tcut=tcut_default,tlow=tlow_default):
        """ Fits the viscosity data to correlation with assigned weights """
        
        eta_data = eta_data[t_data<tcut]
        w8_data = w8_data[t_data<tcut]
        t_data = t_data[t_data<tcut]

        eta_data = eta_data[t_data>tlow]
        w8_data = w8_data[t_data>tlow]
        t_data = t_data[t_data>tlow]
        
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
           
    def gen_GK_MCMC_all(self):
        
        base_path_all = self.base_path_all
        GK_MCMC_all = []
        
        for base_path_i in base_path_all:
            
            GK_MCMC_all.append(GreenKubo_states(base_path=base_path_i))
            
        return GK_MCMC_all
       

class GreenKubo_states():
    def __init__(self,base_path='',fpath_all=None):
        try:
            self.fpath_all = self.gen_fpath_all()
        except:
            self.fpath_all = fpath_all
        self.base_path = base_path
        self.total_path_all = [self.base_path + fpath for fpath in self.fpath_all]    
        self.nStates = len(self.fpath_all)
        self.GK_all = self.gen_GK_all()
        
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
                
    def gen_GK_all(self):
        
        total_path_all = self.total_path_all
        GK_all = []
        
        for path in total_path_all:
            
            GK_all.append(GreenKubo(path))
    
        return GK_all
    
class GreenKubo_SaturatedReplicates():
    def __init__(self,irho,nReps,base_path='',tcut=tcut_default):
        self.base_path = base_path
        self.irho = irho
        self.nReps = nReps
        self.tcut = tcut
        self.fpath_all = self.gen_fpath_all()
        self.total_path_all = [self.base_path + fpath for fpath in self.fpath_all]    
        assert nReps == len(self.fpath_all), 'Number of replicates is incorrect'
        self.GK_Reps = self.gen_GK_Reps()
        self.nTime = len(self.GK_Reps[0].t_GK)    
        self.GK_Reps_all, self.t_GK = self.compile_GK_Reps() 
        self.w8_model = self.w8_hat()                     
        self.GK_Reps_avg = self.GK_time_avg()
        self.eta_inf, self.opt_fit = self.calc_eta_inf()
        self.eta_boots, self.eta_low, self.eta_high, self.opt_fit_boots, self.tcut_boots, self.b_boots = self.bootstrap_eta()
        
    def gen_fpath_all(self):
        #Read in the simulation specifications
    
        irho, nReps = self.irho, self.nReps
    
        fpath_all = []
           
        for iRep in range(nReps):
            
            fpath_all.append('Saturated/rho'+str(irho)+'/Rep'+str(iRep)+'/NVT_eq/NVT_prod/NVT_vis/')
    
        return fpath_all
                
    def gen_GK_Reps(self):
        
        total_path_all = self.total_path_all
        GK_Reps = []
        
        for path in total_path_all:
            
            GK_Reps.append(GreenKubo(path))
    
        return GK_Reps
        
    def w8_hat(self,fit=False):
        t_GK = self.t_GK
        if not fit:
            w8_model = t_GK**(-0.5)
        return w8_model
            
    def compile_GK_Reps(self):
        """ Compile the Green-Kubo values for all replicates """
        nReps,nTime = self.nReps, self.nTime
        
        GK_Reps = self.GK_Reps
        
        GK_Reps_all = np.zeros([nTime,nReps])
        
        for iRep in range(nReps):
            
            GK_Reps_all[:,iRep] = GK_Reps[iRep].visc_GK
                                
        t_GK = GK_Reps[0].t_GK
            
        return GK_Reps_all, t_GK
    
    def GK_time_avg(self):
        """ Averages the viscosity at each time for different replicates """
        GK_Reps_all = self.GK_Reps_all
        
        GK_Reps_avg = np.mean(GK_Reps_all,axis=1)

        np.savetxt('GK_Reps_avg',GK_Reps_avg.T,fmt='%0.7f')
        
        return GK_Reps_avg            
    
    def calc_eta_inf(self):
        
        GK_Reps_avg,t_GK, w8_model = self.GK_Reps_avg,self.t_GK,self.w8_model
      
        opt_fit = self.fit_eta(t_GK,GK_Reps_avg,w8_model)
        eta_inf = self.calc_eta_estimate(opt_fit)
        
#        np.savetxt('GK_eta_inf',eta_inf,fmt='%0.7f')
        
        return eta_inf, opt_fit
    
    def eta_hat(self,t,params):
        """ Model suggested by Maginn for fitting viscosity """
        
        A = params[0]
        alpha = params[1]
        tau1 = params[2]
        tau2 = params[3]
        
        eta_t = A*alpha*tau1*(1-np.exp(-t/tau1))+A*(1-alpha)*tau2*(1-np.exp(-t/tau2))
        return eta_t
    
    def fit_eta(self,t_data,eta_data,w8_data,tcut=tcut_default,tlow=tlow_default):
        """ Fits the viscosity data to correlation with assigned weights """
        
        eta_data = eta_data[t_data<tcut]
        w8_data = w8_data[t_data<tcut]
        t_data = t_data[t_data<tcut]

        eta_data = eta_data[t_data>tlow]
        w8_data = w8_data[t_data>tlow]
        t_data = t_data[t_data>tlow]
        
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
    
    def bootstrap_eta(self):
        
        GK_Reps_avg,t_GK, w8_model,GK_Reps_all = self.GK_Reps_avg,self.t_GK,self.w8_model, self.GK_Reps_all
        
        nBoots = 1000
        self.nBoots = nBoots
        eta_boots = np.zeros(nBoots)
        tcut_boots = np.zeros(nBoots)
        b_boots = np.zeros(nBoots)
        opt_fit_boots = []
                
        tcut_range = np.linspace(0.8*self.tcut,10*self.tcut)
        b_low = 0.4
        b_high = 0.6
        
        for iBoots in range(nBoots):
            
#                w8_boots = np.random.random(self.nReps) #This has a random weight for the different runs, essentially it randomly samples which runs to include 
#                w8_boots = np.random.randint(0,2,self.nReps) # The problem with this method is that you might have all 0's
            w8_boots = np.random.randint(0,self.nReps,self.nReps) # This method tries to ensure that a large number of runs are included
            while np.sum(w8_boots) == 0: w8_boots = np.random.randint(0,2,self.nReps)
 
            eta_avg = np.sum(GK_Reps_all*w8_boots,axis=1)/np.sum(w8_boots)
#            eta_avg = GK_Reps_avg # This approach does not do any random sampling of replicates
            
            tcut = np.random.choice(tcut_range)
            b_random = np.random.uniform(b_low,b_high)
            w8_model = t_GK**(-b_random)
        
            opt_fit = self.fit_eta(t_GK,eta_avg,w8_model,tcut)
            eta_boots[iBoots] = self.calc_eta_estimate(opt_fit)
            
            opt_fit_boots.append(opt_fit)
            tcut_boots[iBoots] = tcut
            b_boots[iBoots] = b_random
        
        eta_boots_sorted = np.sort(np.array(eta_boots))
        
        ilow = int(nBoots*2.5/100.)
        ihigh = int(nBoots*97.5/100.)
        
        eta_low = eta_boots_sorted[ilow]
        eta_high = eta_boots_sorted[ihigh]
            
        np.savetxt('GK_eta_boots',eta_boots,fmt='%0.7f')
        np.savetxt('GK_tcut_boots',tcut_boots,fmt='%0.7f')
        np.savetxt('GK_b_boots',b_boots,fmt='%0.7f')
            
        return eta_boots, eta_low, eta_high, opt_fit_boots, tcut_boots, b_boots
    
    def plot_bootstraps(self):
        """ Plots the bootstrapped histogram """
        eta_boots, eta_low, eta_high, eta_inf, tcut_boots, b_boots = self.eta_boots, self.eta_low, self.eta_high, self.eta_inf, self.tcut_boots, self.b_boots
        
        nbins = int(len(eta_boots.T)/5)
            
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
        plt.plot([eta_inf,eta_inf],[0,hist_max],'r-',label='Estimated Viscosity')
        plt.xlabel('Viscosity (cP)')
        plt.ylabel('Probability Density (1/cP)')
        plt.legend()
        
        fig.savefig('GK_bootstraps_Reps.pdf') 
        plt.close()
               
#        fig, ax1 = plt.subplots(figsize=(6,6))
#        ax1.plot(tcut_boots, eta_boots, 'b-')
#        ax1.set_ylabel('Viscosity (cP)')
#        ax1.set_xlabel('Cut-off Time (ps)', color='b')
#        ax1.tick_params('x', colors='b')
#        
#        ax2 = ax1.twiny()
#        ax2.plot(b_boots, eta_boots, 'r--')
#        ax2.set_xlabel('Weight Exponent', color='r')
#        ax2.tick_params('x', colors='r')
#        
#        fig.tight_layout()        
#        fig.savefig('GK_tcut_b_boots_Reps.pdf') 
#        plt.close()
        
#        fig = plt.figure(figsize=(6,6))
#        
#        p = plt.scatter(tcut_boots,b_boots,c=eta_boots, cmap='rainbow')
#        cb = plt.colorbar(p)
#        cb.set_label('Viscosity (cP)')
#        plt.xlabel('Cut-off Time (ps)')
#        plt.ylabel('Weight Exponent')
#        
#        fig.savefig('GK_heatmap_boots.pdf')
#        plt.close()
        

    def plot_GK_Reps(self):
        """ Plots all the average GK on the same plot """
        
        nReps, GK_Reps_all,GK_Reps_avg, t_GK, tcut = self.nReps,self.GK_Reps_all,self.GK_Reps_avg, self.t_GK, self.tcut
        
        tplot = np.linspace(0,t_GK.max(),10000)
        GK_plot = self.eta_hat(tplot,self.opt_fit)
        
        GK_boot_low = np.ones(len(tplot))*100.
        GK_boot_high = np.ones(len(tplot))*(-100.)
        
        fig = plt.figure(figsize=(6,6))
        
        color_scheme = ['b','r','g','c','m','y','brown','orange','purple','pink','teal','turquoise','violet','indigo','grey','navy','darkgreen','goldenrod','lawngreen','royalblue','tomato', 'maroon']

#        for iRep in range(nReps):
#            
#            GK_data = GK_Reps_all[:,iRep]
#            
#            plt.plot(t_GK,GK_data,':',alpha=0.4)
        
        plt.plot(t_GK,GK_Reps_all,':',alpha=0.4)
        plt.plot(t_GK,GK_Reps_avg,'k-',label='Average')
        plt.plot(tplot,GK_plot,'b--',label='Fit')
        
        for iBoot in range(self.nBoots):
            
            GK_boot = self.eta_hat(tplot,self.opt_fit_boots[iBoot])
            
            for it in range(len(tplot)):
                
                if GK_boot[it] > GK_boot_high[it]:
                    GK_boot_high[it] = GK_boot[it]
                    
                if GK_boot[it] < GK_boot_low[it]:
                    GK_boot_low[it] = GK_boot[it]
                
#            plt.plot(tplot,GK_boot,'r-',alpha=0.2)
        
        plt.plot(tplot,GK_boot_low,'r-.')
        plt.plot(tplot,GK_boot_high,'r-.',label='Bootstraps')
        plt.plot([tcut,tcut],[0,GK_plot.max()],'g-.',label='Cut-off')
        
        plt.xlabel('Time (ps)')
#        plt.xlim([0,tcut])
        
        plt.ylim(ymin=0)
        plt.ylabel('Viscosity (cP)')
        plt.legend()
                
        fig.savefig('GK_Reps_all.pdf') 
        plt.close()    
        
    def plot_scan_tcut_b(self):
        
        GK_Reps_avg,t_GK = self.GK_Reps_avg,self.t_GK
        
        eta_avg = GK_Reps_avg # This approach does not do any random sampling of replicates
                
        tcut_range = np.linspace(0.8*self.tcut,10*self.tcut,50)
        b_range = np.linspace(0.4,0.6,50)
        
        eta_scan = np.zeros([len(tcut_range),len(b_range)])
        
#        for it, tcut_i in enumerate(tcut_range):
#            
#            for ib, b_i in enumerate(b_range):
#            
#                w8_model = t_GK**(-b_i)
#            
#                opt_fit = self.fit_eta(t_GK,eta_avg,w8_model,tcut_i)
#                eta_scan[it,ib] = self.calc_eta_estimate(opt_fit)
#                        
#        np.savetxt('GK_eta_scan',eta_scan,fmt='%0.7f')
#        
#        f = plt.figure()
#        plt.contour(tcut_range,b_range,eta_scan)
#        plt.ylabel('Weight Exponent')
#        plt.xlabel('Cut-off Time (ps)')
#        plt.ylim([min(b_range),max(b_range)])
#        plt.xlim([min(tcut_range),max(tcut_range)])
#        plt.colorbar()
#        f.savefig('eta_scan_tcut_b.pdf')
        
        tcut_range = np.linspace(8,500,200)
        
        w8_model = t_GK**(-0.5)
                
        eta_scan_tcut = np.zeros(len(tcut_range))
        
        for it, tcut_i in enumerate(tcut_range):

            opt_fit = self.fit_eta(t_GK,eta_avg,w8_model,tcut_i)
            eta_scan_tcut[it] = self.calc_eta_estimate(opt_fit)
                    
        np.savetxt('GK_eta_scan_tcut',eta_scan_tcut,fmt='%0.7f')
        
        f = plt.figure()
        plt.plot(tcut_range,eta_scan_tcut)
        plt.ylabel('Viscosity (cP)')
        plt.xlabel('Cut-off Time (ps)')
        plt.xlim([min(tcut_range),max(tcut_range)])
        f.savefig('eta_scan_tcut.pdf')
        
        b_range = np.linspace(0.2,0.8,200)
        
        eta_scan_b = np.zeros(len(tcut_range))
        
        for ib, b_i in enumerate(b_range):
            
            w8_model = t_GK**(-b_i)

            opt_fit = self.fit_eta(t_GK,eta_avg,w8_model,500.)
            eta_scan_b[ib] = self.calc_eta_estimate(opt_fit)
                    
        np.savetxt('GK_eta_scan_b',eta_scan_b,fmt='%0.7f')
        
        f = plt.figure()
        plt.plot(b_range,eta_scan_b)
        plt.ylabel('Viscosity (cP)')
        plt.xlabel('Weight Exponent')
        plt.xlim([min(b_range),max(b_range)])
        f.savefig('eta_scan_b.pdf')
            
        return eta_scan, tcut_range, b_range
    
class GreenKubo_SaturatedMCMC():
    def __init__(self,irho,ilow,ihigh,nReps=1,tcut=tcut_default):
        
        self.ilow, self.ihigh = ilow, ihigh
        self.irho = irho
        self.nMCMC = (ihigh-ilow+1)
        self.nReps = nReps
        self.tcut = tcut
        self.fpath_all = self.gen_fpath_all()  
        assert self.nMCMC*self.nReps == len(self.fpath_all), 'Number of MCMC/replicates is incorrect'
        self.GK_MCMC = self.gen_GK_MCMC()
        self.nTime = len(self.GK_MCMC[0].t_GK)    
        self.GK_MCMC_all, self.t_GK = self.compile_GK_MCMC() 
        self.w8_model = self.w8_hat()                     
        self.GK_MCMC_avg = self.GK_time_avg()
        self.eta_inf, self.opt_fit = self.calc_eta_inf()
#        self.eta_boots, self.eta_low, self.eta_high, self.opt_fit_boots, self.tcut_boots, self.b_boots = self.bootstrap_eta()
        self.eta_boots, self.eta_low, self.eta_high, self.opt_fit_boots = self.bootstrap_eta_alt()
        
    def gen_fpath_all(self):
        #Read in the simulation specifications
    
        irho, nReps,nMCMC = self.irho, self.nReps,self.nMCMC
    
        fpath_all = []
           
        for iMCMC in range(nMCMC):
            
            for iRep in range(nReps):
            
                fpath_all.append('MCMC_'+str(iMCMC)+'/Saturated/rho'+str(irho)+'/Rep'+str(iRep)+'/NVT_eq/NVT_prod/NVT_vis/')
    
        return fpath_all
                
    def gen_GK_MCMC(self):
        
        fpath_all = self.fpath_all
        GK_MCMC = []
        
        for path in fpath_all:
            
            GK_MCMC.append(GreenKubo(path))
    
        return GK_MCMC
        
    def w8_hat(self,fit=False):
        t_GK = self.t_GK
        if not fit:
            w8_model = t_GK**(-0.5)
        return w8_model
            
    def compile_GK_MCMC(self):
        """ Compile the Green-Kubo values for all MCMC samples """
        nMCMC,nTime,GK_MCMC = self.nMCMC, self.nTime, self.GK_MCMC
        
        GK_MCMC_all = np.zeros([nTime,nMCMC])
        # Currently this does not allow for multiple replicates. If I end up needing that I can add that feature
        for iMCMC in range(nMCMC):
            
            GK_MCMC_all[:,iMCMC] = GK_MCMC[iMCMC].visc_GK
                                
        t_GK = GK_MCMC[0].t_GK
            
        return GK_MCMC_all, t_GK
    
    def GK_time_avg(self):
        """ Averages the viscosity at each time for different MCMC samples """
        GK_MCMC_all, irho = self.GK_MCMC_all, self.irho
        
        GK_MCMC_avg = np.mean(GK_MCMC_all,axis=1)

        np.savetxt('GK_MCMC_avg_rho'+str(irho),GK_MCMC_avg.T,fmt='%0.7f')
        
        return GK_MCMC_avg            
    
    def calc_eta_inf(self):
        
        GK_MCMC_avg,t_GK, w8_model = self.GK_MCMC_avg,self.t_GK,self.w8_model
      
        opt_fit = self.fit_eta(t_GK,GK_MCMC_avg,w8_model)
        eta_inf = self.calc_eta_estimate(opt_fit)
        
#        np.savetxt('GK_eta_inf',eta_inf,fmt='%0.7f')
        
        return eta_inf, opt_fit
    
    def eta_hat(self,t,params):
        """ Model suggested by Maginn for fitting viscosity """
        
        A = params[0]
        alpha = params[1]
        tau1 = params[2]
        tau2 = params[3]
        
        eta_t = A*alpha*tau1*(1-np.exp(-t/tau1))+A*(1-alpha)*tau2*(1-np.exp(-t/tau2))
        return eta_t
    
    def fit_eta(self,t_data,eta_data,w8_data,tcut=tcut_default,tlow=tlow_default):
        """ Fits the viscosity data to correlation with assigned weights """
        
        eta_data = eta_data[t_data<tcut]
        w8_data = w8_data[t_data<tcut]
        t_data = t_data[t_data<tcut]

        eta_data = eta_data[t_data>tlow]
        w8_data = w8_data[t_data>tlow]
        t_data = t_data[t_data>tlow]
        
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
    
    def bootstrap_eta(self):
        
        GK_MCMC_avg,t_GK, w8_model,GK_MCMC_all = self.GK_MCMC_avg,self.t_GK,self.w8_model, self.GK_MCMC_all
        
        nBoots = 1000
        self.nBoots = nBoots
        eta_boots = np.zeros(nBoots)
        tcut_boots = np.zeros(nBoots)
        b_boots = np.zeros(nBoots)
        opt_fit_boots = []
                
        tcut_range = np.linspace(0.8*self.tcut,10*self.tcut)
        b_low = 0.4
        b_high = 0.6
        
        for iBoots in range(nBoots):
            
#                w8_boots = np.random.random(self.nMCMC) #This has a random weight for the different runs, essentially it randomly samples which runs to include 
#                w8_boots = np.random.randint(0,2,self.nMCMC) # The problem with this method is that you might have all 0's
            w8_boots = np.random.randint(0,self.nMCMC,self.nMCMC) # This method tries to ensure that a large number of runs are included
            while np.sum(w8_boots) == 0: w8_boots = np.random.randint(0,2,self.nMCMC)
 
            eta_avg = np.sum(GK_MCMC_all*w8_boots,axis=1)/np.sum(w8_boots)
#            eta_avg = GK_MCMC_avg # This approach does not do any random sampling of replicates
            
            tcut = np.random.choice(tcut_range)
            b_random = np.random.uniform(b_low,b_high)
            w8_model = t_GK**(-b_random)
        
            opt_fit = self.fit_eta(t_GK,eta_avg,w8_model,tcut)
            eta_boots[iBoots] = self.calc_eta_estimate(opt_fit)
            
            opt_fit_boots.append(opt_fit)
            tcut_boots[iBoots] = tcut
            b_boots[iBoots] = b_random
        
        eta_boots_sorted = np.sort(np.array(eta_boots))
        
        ilow = int(nBoots*2.5/100.)
        ihigh = int(nBoots*97.5/100.)
        
        eta_low = eta_boots_sorted[ilow]
        eta_high = eta_boots_sorted[ihigh]
            
        np.savetxt('GK_eta_boots',eta_boots,fmt='%0.7f')
        np.savetxt('GK_tcut_boots',tcut_boots,fmt='%0.7f')
        np.savetxt('GK_b_boots',b_boots,fmt='%0.7f')
            
        return eta_boots, eta_low, eta_high, opt_fit_boots, tcut_boots, b_boots

    def bootstrap_eta_alt(self):
        
        GK_MCMC_avg,t_GK, w8_model,GK_MCMC_all,irho = self.GK_MCMC_avg,self.t_GK,self.w8_model, self.GK_MCMC_all, self.irho
        
        nBootsMCMC = 100
        nBootsReps = 200
        
        nBoots = 200
        self.nBoots = nBoots
        eta_boots = np.zeros(nBoots*nBootsMCMC)
        opt_fit_boots = []
                
        tcut_range = np.linspace(0.8*self.tcut,10*self.tcut)
        b_low = 0.4
        b_high = 0.6

        ieta = 0
        
        for iBootsMCMC in range(nBootsMCMC):
                       
           GK_MCMC_reps = np.zeros([GK_MCMC_all.shape[0],nBootsReps])
           print('iBootsMCMC = '+str(iBootsMCMC))
           for iBootsReps in range(nBootsReps):
               
               GK_MCMC_reps[:,iBootsReps] = GK_MCMC_all[:,np.random.randint(0,GK_MCMC_all.shape[1])]
               
           for iBoots in range(nBoots):
                
    #                w8_boots = np.random.random(self.nMCMC) #This has a random weight for the different runs, essentially it randomly samples which runs to include 
    #                w8_boots = np.random.randint(0,2,self.nMCMC) # The problem with this method is that you might have all 0's
                w8_boots = np.random.randint(0,nBootsReps,nBootsReps) # This method tries to ensure that a large number of runs are included
                while np.sum(w8_boots) == 0: w8_boots = np.random.randint(0,2,nBootsReps)
     
                eta_avg = np.sum(GK_MCMC_reps*w8_boots,axis=1)/np.sum(w8_boots)
    #            eta_avg = GK_MCMC_avg # This approach does not do any random sampling of replicates
                
                tcut = np.random.choice(tcut_range)
                b_random = np.random.uniform(b_low,b_high)
                w8_model = t_GK**(-b_random)
            
                opt_fit = self.fit_eta(t_GK,eta_avg,w8_model,tcut)
                eta_boots[ieta] = self.calc_eta_estimate(opt_fit)
                
                opt_fit_boots.append(opt_fit)
                ieta += 1

        eta_boots_sorted = np.sort(np.array(eta_boots))
        
        ilow = int(len(eta_boots)*2.5/100.)
        ihigh = int(len(eta_boots)*97.5/100.)
        
        eta_low = eta_boots_sorted[ilow]
        eta_high = eta_boots_sorted[ihigh]
        
        np.savetxt('GK_eta_boots_rho'+str(irho),eta_boots,fmt='%0.7f')

        return eta_boots, eta_low, eta_high, opt_fit_boots
    
    def plot_bootstraps(self):
        """ Plots the bootstrapped histogram """
        eta_boots, eta_low, eta_high, eta_inf, irho = self.eta_boots, self.eta_low, self.eta_high, self.eta_inf, self.irho
        
        try: 
            tcut_boots, b_boots = self.tcut_boots, self.b_boots 
        except: 
            pass
        
        nbins = int(len(eta_boots.T)/5)
            
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
        plt.plot([eta_inf,eta_inf],[0,hist_max],'r-',label='Estimated Viscosity')
        plt.xlabel('Viscosity (cP)')
        plt.ylabel('Probability Density (1/cP)')
        plt.legend()
        
        fig.savefig('GK_bootstraps_MCMC_rho'+str(irho)+'.pdf') 
        plt.close()
               
#        fig, ax1 = plt.subplots(figsize=(6,6))
#        ax1.plot(tcut_boots, eta_boots, 'b-')
#        ax1.set_ylabel('Viscosity (cP)')
#        ax1.set_xlabel('Cut-off Time (ps)', color='b')
#        ax1.tick_params('x', colors='b')
#        
#        ax2 = ax1.twiny()
#        ax2.plot(b_boots, eta_boots, 'r--')
#        ax2.set_xlabel('Weight Exponent', color='r')
#        ax2.tick_params('x', colors='r')
#        
#        fig.tight_layout()        
#        fig.savefig('GK_tcut_b_boots_MCMC.pdf') 
#        plt.close()
        
#        fig = plt.figure(figsize=(6,6))
#        
#        p = plt.scatter(tcut_boots,b_boots,c=eta_boots, cmap='rainbow')
#        cb = plt.colorbar(p)
#        cb.set_label('Viscosity (cP)')
#        plt.xlabel('Cut-off Time (ps)')
#        plt.ylabel('Weight Exponent')
#        
#        fig.savefig('GK_heatmap_boots.pdf')
#        plt.close()
        

    def plot_GK_MCMC(self):
        """ Plots all the average GK on the same plot """
        
        nMCMC, GK_MCMC_all,GK_MCMC_avg, t_GK, tcut, irho = self.nMCMC,self.GK_MCMC_all,self.GK_MCMC_avg, self.t_GK, self.tcut, self.irho
        
        tplot = np.linspace(0,t_GK.max(),10000)
        GK_plot = self.eta_hat(tplot,self.opt_fit)
        
        GK_boot_low = np.ones(len(tplot))*100.
        GK_boot_high = np.ones(len(tplot))*(-100.)
        
        fig = plt.figure(figsize=(6,6))
        
        color_scheme = ['b','r','g','c','m','y','brown','orange','purple','pink','teal','turquoise','violet','indigo','grey','navy','darkgreen','goldenrod','lawngreen','royalblue','tomato', 'maroon']

#        for iMCMC in range(nMCMC):
#            
#            GK_data = GK_MCMC_all[:,iMCMC]
#            
#            plt.plot(t_GK,GK_data,':',alpha=0.4)
        
        plt.plot(t_GK,GK_MCMC_all,':',alpha=0.4)
        plt.plot(t_GK,GK_MCMC_avg,'k-',label='Average')
        plt.plot(tplot,GK_plot,'b--',label='Fit')
        
        for iBoot in range(self.nBoots):
            
            GK_boot = self.eta_hat(tplot,self.opt_fit_boots[iBoot])
            
            for it in range(len(tplot)):
                
                if GK_boot[it] > GK_boot_high[it]:
                    GK_boot_high[it] = GK_boot[it]
                    
                if GK_boot[it] < GK_boot_low[it]:
                    GK_boot_low[it] = GK_boot[it]
                
#            plt.plot(tplot,GK_boot,'r-',alpha=0.2)
        
        plt.plot(tplot,GK_boot_low,'r-.')
        plt.plot(tplot,GK_boot_high,'r-.',label='Bootstraps')
        plt.plot([tcut,tcut],[0,GK_plot.max()],'g-.',label='Cut-off')
        
        plt.xlabel('Time (ps)')
#        plt.xlim([0,tcut])
        
        plt.ylim(ymin=0)
        plt.ylabel('Viscosity (cP)')
        plt.legend()
                
        fig.savefig('GK_MCMC_all_rho'+str(irho)+'.png') 
        plt.close()    
        
    def plot_scan_tcut_b(self):
        
        GK_MCMC_avg,t_GK = self.GK_MCMC_avg,self.t_GK
        
        eta_avg = GK_MCMC_avg # This approach does not do any random sampling of replicates
                
        tcut_range = np.linspace(0.8*self.tcut,10*self.tcut,50)
        b_range = np.linspace(0.4,0.6,50)
        
        eta_scan = np.zeros([len(tcut_range),len(b_range)])
        
#        for it, tcut_i in enumerate(tcut_range):
#            
#            for ib, b_i in enumerate(b_range):
#            
#                w8_model = t_GK**(-b_i)
#            
#                opt_fit = self.fit_eta(t_GK,eta_avg,w8_model,tcut_i)
#                eta_scan[it,ib] = self.calc_eta_estimate(opt_fit)
#                        
#        np.savetxt('GK_eta_scan',eta_scan,fmt='%0.7f')
#        
#        f = plt.figure()
#        plt.contour(tcut_range,b_range,eta_scan)
#        plt.ylabel('Weight Exponent')
#        plt.xlabel('Cut-off Time (ps)')
#        plt.ylim([min(b_range),max(b_range)])
#        plt.xlim([min(tcut_range),max(tcut_range)])
#        plt.colorbar()
#        f.savefig('eta_scan_tcut_b.pdf')
        
        tcut_range = np.linspace(8,500,200)
        
        w8_model = t_GK**(-0.5)
                
        eta_scan_tcut = np.zeros(len(tcut_range))
        
        for it, tcut_i in enumerate(tcut_range):

            opt_fit = self.fit_eta(t_GK,eta_avg,w8_model,tcut_i)
            eta_scan_tcut[it] = self.calc_eta_estimate(opt_fit)
                    
        np.savetxt('GK_eta_scan_tcut',eta_scan_tcut,fmt='%0.7f')
        
        f = plt.figure()
        plt.plot(tcut_range,eta_scan_tcut)
        plt.ylabel('Viscosity (cP)')
        plt.xlabel('Cut-off Time (ps)')
        plt.xlim([min(tcut_range),max(tcut_range)])
        f.savefig('eta_scan_tcut.pdf')
        
        b_range = np.linspace(0.2,0.8,200)
        
        eta_scan_b = np.zeros(len(tcut_range))
        
        for ib, b_i in enumerate(b_range):
            
            w8_model = t_GK**(-b_i)

            opt_fit = self.fit_eta(t_GK,eta_avg,w8_model,500.)
            eta_scan_b[ib] = self.calc_eta_estimate(opt_fit)
                    
        np.savetxt('GK_eta_scan_b',eta_scan_b,fmt='%0.7f')
        
        f = plt.figure()
        plt.plot(b_range,eta_scan_b)
        plt.ylabel('Viscosity (cP)')
        plt.xlabel('Weight Exponent')
        plt.xlim([min(b_range),max(b_range)])
        f.savefig('eta_scan_b.pdf')
            
        return eta_scan, tcut_range, b_range

class GreenKubo():
    def __init__(self, fpath):
        self.fpath = fpath
        self.t_GK, self.visc_GK = self.load_visc_GK()
        
    def load_visc_GK(self):
        """ Reads in the viscosities for the Green-Kubo """
        fpath = self.fpath
        g_start = 25
        t_visc_GK = open(fpath+'visco.xvg','r').readlines()[g_start:] #Read all lines starting at g_start
        n_t = len(t_visc_GK)
                    
        t_GK = np.zeros(n_t)
        visc_GK = np.zeros(n_t)
         
        for i_t in range(n_t):
            t_GK[i_t] = float(t_visc_GK[i_t].split()[0])
            visc_GK[i_t] = float(t_visc_GK[i_t].split()[1])
        return t_GK, visc_GK
    
def main():
       
    parser = argparse.ArgumentParser()
    parser.add_argument("-il","--ilow",type=int,help="Specify the lowest iMCMC value")
    parser.add_argument("-ih","--ihigh",type=int,help="Specify the high iMCMC value")
    parser.add_argument("-nR","--nReps",type=int,help="Specify the number of replicates")
    parser.add_argument("--sat",action='store_true',help="If just a single saturated condition")
    parser.add_argument("--scan",action='store_true',help="If you want to scan tcut and b parameter space")
    parser.add_argument("-ir","--irho",type=int,help="Specify the irho")
    args = parser.parse_args()
    
    if args.sat:
        
        if args.ihigh:
            
            MCMC_Mie = GreenKubo_SaturatedMCMC(args.irho,args.ilow,args.ihigh,args.nReps)
            MCMC_Mie.plot_GK_MCMC()
            
        elif args.nReps:
        
            MCMC_Mie = GreenKubo_SaturatedReplicates(args.irho,args.nReps)
            MCMC_Mie.plot_GK_Reps()            
        
        MCMC_Mie.plot_bootstraps()
        
        if args.scan:
        
            MCMC_Mie.plot_scan_tcut_b()
        
    else:

        MCMC_Mie = GreenKubo_MCMC(args.ilow,args.ihigh)
        MCMC_Mie.plot_GK_all()
        MCMC_Mie.plot_eta_sat()
        MCMC_Mie.plot_eta_IT()
        MCMC_Mie.plot_bootstraps()
   
    
if __name__ == '__main__':
    '''
    python GreenKubo_analyze.py --ilow int --ihigh int --nReps int --irho int --sat --scan
  
    '''

    main()   