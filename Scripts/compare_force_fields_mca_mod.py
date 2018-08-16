"""
Compare simulation results with experimental data in TDE and REFPROP correlations
for several force fields and compounds.
Must call from ~ with the following files already created:
    SaturatedSettings.txt
    REFPROP_eta_sat.txt
    TDE_eta_sat.txt
    GK_eta_boots_rho"$irho"
Must run from ~\ home
Example: python compare_force_fields.py --comp C3H8 C4H10 C8H18 --mod Potoff TraPPE TAMie
"""

from __future__ import division
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import os

font = {'size' : '24'}
plt.rc('font',**font)

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--comp",type=str,nargs='+',choices=['C3H8','C4H10','C8H18','C12H26','IC4H10','IC5H12','IC6H14','IC8H18','NEOC5H12', '23DMButane', '3MPentane', 'C16H34','C22H46'],help="Specify the compound(s) to analyze")
    parser.add_argument("--devPlot",action='store_true',help="If deviation plots are desired")
    parser.add_argument("--uncer",action='store_true',help="If REFPROP uncertainty should be plotted")
    parser.add_argument("--TDEuncer",action='store_true',help="If TDE uncertainty should be plotted")
    parser.add_argument("-m","--mod",type=str,nargs='+',choices=['Potoff','TraPPE','TAMie','AUA4','Exp6','Mie_14_Bayesian_Viscosity_200_MCMC','Mie_16_Bayesian_Viscosity_200_MCMC'],help="Specify the model(s) to analyze")
    parser.add_argument("-v","--validation",type=str,nargs='+',choices=['TAMie','TraPPE','Potoff'],help="Validation statepoints to plot using any force field.")
    parser.add_argument("--settings",type=int,nargs='?',help="Settings for the location of labels and the REFPROP uncertainty.")
    args = parser.parse_args()
    
    if args.comp and args.mod:
        
        model_list=args.mod
        compound_list=args.comp
        
        #color_list = {'C3H8':'r','C4H10':'b','C8H18':'g','IC8H18':'m'}
        color_list = {'Potoff':'r','TraPPE':'b','TAMie':'g','AUA4':'m','Mie_14_Bayesian_Viscosity_200_MCMC':'b','Mie_16_Bayesian_Viscosity_200_MCMC':'r'}
        line_list = {'C3H8':'-','C4H10':'--','C8H18':'-.','IC8H18': ':'}
        symbol_list = {'Potoff':'o','TraPPE':'s','TAMie':'^','AUA4':'v','Mie_14_Bayesian_Viscosity_200_MCMC':'s','Mie_16_Bayesian_Viscosity_200_MCMC':'o'}
        tail_path_list = {'C3H8':'_N400','C4H10':'_N400','C8H18':'_N400','C12H26':'_N400','IC8H18':'_N400', 'C16H34':'_N200','3MPentane':'_N400', 'NEOC5H12':'_N400','IC6H14':'_N400','23DMButane':'_N400','IC5H12':'_N400','IC4H10':'_N400','C22H46':'_N200'}
        label_list = {'C3H8':r'$n$-C$_3$','C4H10':r'$n$-C$_4$','C8H18':r'$n$-C$_8$','C12H26':r'$n$-C$_{12}$','IC4H10':r'$i$-C$_4$','IC5H12':r'$i$-C$_5$','IC6H14':r'$i$-C$_6$','IC8H18':r'$i$-C$_8$','NEOC5H12':r'neo-C$_5$','3MPentane':r'3-MC$_5$','23DMButane':r'2,3-DMC$_4$', 'C16H34':r'$n$-C$_1$$_6$','C22H46':r'$n$-C$_2$$_2$'}
        ### For cP use first list, for Pa-s use second list
        #label_position_list = {'C3H8':[330,0.03],'C4H10':[400,0.03],'C8H18':[530,0.03],'IC4H10':[380,0.035],'IC5H12':[420,0.03],'IC6H14':[470,0.03],'IC8H18':[520,0.03],'NEOC5H12':[390,0.1]}
        label_position_list = {'3MPentane':[250,.00003],'23DMButane':[310,.00003],'C3H8':[330,0.00003],'C4H10':[400,0.00003],'C8H18':[530,0.00003],'C12H26':[630,0.00003],'IC4H10':[380,0.000035],'IC5H12':[420,0.00003],'IC6H14':[470,0.00003],'IC8H18':[520,0.00003],'NEOC5H12':[390,0.0001],'C16H34':[650,0.0002],'C22H46':[690,0.0002]}
        nrhomax = {'C3H8':10,'C4H10':5,'C8H18':5,'C12H26':5,'IC4H10':5,'IC5H12':5,'IC6H14':5,'IC8H18':5,'NEOC5H12':5,'3MPentane':5,'23DMButane':5,'C16H34':5,'C22H46':5}
        REF_uncer = {'C3H8':.005,'C4H10':.005,'C8H18':.005,'C12H26':.005,'IC4H10':.03,'IC5H12':.04,'IC6H14':.02,'IC8H18':.05,'NEOC5H12':.10,'3MPentane':.02,'23DMButane':.05,'C16H34':.02,'C22H46':.10}
    
        # Machinery for implementing vertical shifts of various compounds. 
        shifts = {'C3H8':0.0,'C4H10':0.0,'C8H18':0.0,'C12H26':0.0,'IC4H10':0.0,'IC5H12':0.0,'IC6H14':0.0,'IC8H18':0.0,'NEOC5H12':0.0,'3MPentane':0.0,'23DMButane':0.0,'C16H34':0.0,'C22H46':0.0}
        # Short branched alkanes plot iC4 iC5 23DMC4 neoC5
        if args.settings == 1:
            shifts = {'C3H8':0.0,'C4H10':1.0,'C8H18':1.0,'C12H26':1.0,'IC4H10':0.0,'IC5H12':0.5,'IC6H14':1.0,'IC8H18':1.0,'NEOC5H12':-.5,'3MPentane':-1.0,'23DMButane':1.0,'C16H34':1.0,'C22H46':1.0}
            label_position_list = {'3MPentane':[250,.00003],'23DMButane':[430,-2.95],'C3H8':[330,0.00003],'C4H10':[400,0.00003],'C8H18':[530,0.00003],'C12H26':[630,0.00003],'IC4H10':[410,-4.4],'IC5H12':[470,-3.95],'IC6H14':[470,0.00003],'IC8H18':[520,0.00003],'NEOC5H12':[435,-4.95],'C16H34':[650,0.0002],'C22H46':[690,0.0002]}
        # Long: 3MPentane IC8 IC6
        elif args.settings == 2:
            shifts = {'C3H8':0.0,'C4H10':1.0,'C8H18':1.0,'C12H26':1.0,'IC4H10':0.0,'IC5H12':0.5,'IC6H14':0.0,'IC8H18':1.0,'NEOC5H12':-.5,'3MPentane':-1.0,'23DMButane':1.0,'C16H34':1.0,'C22H46':1.0}
            label_position_list = {'3MPentane':[500,-5.3],'23DMButane':[430,-3.6],'C3H8':[330,0.00003],'C4H10':[400,0.00003],'C8H18':[530,0.00003],'C12H26':[630,0.00003],'IC4H10':[420,-4.5],'IC5H12':[460,-4.0],'IC6H14':[500,-4.3],'IC8H18':[505,-3.0],'NEOC5H12':[435,-5.0],'C16H34':[650,0.0002],'C22H46':[690,0.0002]}
        # IC6 3MP 23DMB
        elif args.settings == 3:
            shifts = {'C3H8':0.0,'C4H10':1.0,'C8H18':1.0,'C12H26':1.0,'IC4H10':0.0,'IC5H12':0.5,'IC6H14':-1.0,'IC8H18':1.0,'NEOC5H12':-.5,'3MPentane':0.0,'23DMButane':1.0,'C16H34':1.0,'C22H46':1.0}
            label_position_list = {'3MPentane':[435,-5.0],'23DMButane':[430,-3.7],'C3H8':[330,0.00003],'C4H10':[400,0.00003],'C8H18':[530,0.00003],'C12H26':[630,0.00003],'IC4H10':[420,-4.5],'IC5H12':[460,-4.0],'IC6H14':[460,-4.5],'IC8H18':[505,-3.0],'NEOC5H12':[435,-5.0],'C16H34':[650,0.0002],'C22H46':[690,0.0002]}
        # All normal alkanes, in order from C3 to C8
        elif args.settings == 4:
            shifts = {'C3H8':0.0,'C4H10':0.0,'C8H18':0.0,'C12H26':0.0,'IC4H10':0.0,'IC5H12':0.0,'IC6H14':0.0,'IC8H18':0.0,'NEOC5H12':0.0,'3MPentane':0.0,'23DMButane':0.0,'C16H34':0.0,'C22H46':0.0}
            label_position_list = {'3MPentane':[435,-5.0],'23DMButane':[430,-3.7],'C3H8':[340,-4.5],'C4H10':[420,-4.5],'C8H18':[530,-4.5],'C12H26':[630,0.00003],'IC4H10':[420,-4.5],'IC5H12':[460,-4.0],'IC6H14':[460,-4.5],'IC8H18':[505,-3.0],'NEOC5H12':[435,-5.0],'C16H34':[650,0.0002],'C22H46':[690,0.0002]}
        # All normal alkanes, in order from C12 to C22
        elif args.settings == 5:
            shifts = {'C3H8':0.0,'C4H10':0.0,'C8H18':0.0,'C12H26':0.0,'IC4H10':0.0,'IC5H12':0.0,'IC6H14':0.0,'IC8H18':0.0,'NEOC5H12':-.5,'3MPentane':0.0,'23DMButane':0.0,'C16H34':0.0,'C22H46':0.0}
            label_position_list = {'3MPentane':[435,-5.0],'23DMButane':[430,-3.7],'C3H8':[340,-4.5],'C4H10':[420,-4.5],'C8H18':[530,-4.5],'C12H26':[600,-4.57],'IC4H10':[420,-4.2],'IC5H12':[460,-4.0],'IC6H14':[460,-4.5],'IC8H18':[505,-3.0],'NEOC5H12':[435,-5.0],'C16H34':[670,-4.57],'C22H46':[750,-4.57]}

        

        #eta_avg_all = {0:[],1:[],2:[],3:[],4:[]}
        #eta_95_all = {0:[],1:[],2:[],3:[],4:[]}
        if args.validation:
            vallist = args.validation
        else:
            vallist = []
        
        if args.devPlot:
            fig,axarr = plt.subplots(ncols=1,nrows=2,figsize=[10,20],squeeze=False)
        else:
            fig,axarr = plt.subplots(ncols=1,nrows=1,figsize=[10,10],squeeze=False)
        #axarr[0,0].set_yscale("log")  # UNSET THIS

        axarr[0,0].plot([],[],'kx',markersize=10,markeredgewidth=2,mfc='None',label='TDE Exp. Data')
        axarr[0,0].plot([],[],'k-',linewidth=3,label='REFPROP')

        for model in model_list: axarr[0,0].plot([],[],color_list[model]+symbol_list[model],markersize=10,markeredgewidth=2,mfc='None',label=model)
        
        for iC, compound in enumerate(compound_list):
            # Look for special validation data to plot
            for ic, model in enumerate(vallist):
                root_path = compound+'/Gromacs/'+model+"_Saturation_Viscosity/"
                try:
                    try:
                        sim_sat = np.loadtxt(root_path+model+'_SaturatedSettings.txt',skiprows=1)
                    except:
                        sim_sat = np.loadtxt(root_path+model+'_SaturationSettings.txt',skiprows=1)
                    try:
                        Tsat_sim = sim_sat[:,2]
                    except:
                        Tsat_sim = np.array([sim_sat[2]])
                    RP_Tsat_eta = np.loadtxt(root_path+'REFPROP_eta_sat.txt',skiprows=1)
                    RP_Tsat = RP_Tsat_eta[:,0]
                    RP_eta_sat = RP_Tsat_eta[:,1] / 1000. #[Pa-s]
                    RP_eta_sim = np.interp(Tsat_sim,RP_Tsat,RP_eta_sat)
                except:
                    print("Unable to load validation data for "+compound)
                try:
                    eta_MCMC = np.loadtxt(root_path+'GK_eta_boots_rho0') / 1000.  #[Pa-s]
                    per_error_MCMC = (eta_MCMC-RP_eta_sim[0])/RP_eta_sim[0]*100.
                    eta_MCMC = np.log10(eta_MCMC) + shifts[compound]
                    eta_MCMC_sorted = np.sort(eta_MCMC)
                    i95 = int(0.025*len(eta_MCMC_sorted))
                    eta_MCMC_95 = eta_MCMC_sorted[i95:-i95]
                    eta_MCMC_avg = np.mean(eta_MCMC_95) #Only take average without the outliers          
                    eta_MCMC_low = eta_MCMC_avg - eta_MCMC_95[0]
                    eta_MCMC_high = eta_MCMC_95[-1] - eta_MCMC_avg
                    per_error_MCMC_sorted = np.sort(per_error_MCMC)
                    i95 = int(0.025*len(per_error_MCMC_sorted))
                    per_error_MCMC_plot = per_error_MCMC_sorted[i95:-i95]
                    per_error_MCMC_avg = np.mean(per_error_MCMC_plot) #Only take average without the outliers          
                    per_error_MCMC_low = per_error_MCMC_avg - per_error_MCMC_plot[0]
                    per_error_MCMC_high = per_error_MCMC_plot[-1] - per_error_MCMC_avg

                    print(Tsat_sim[0])
                    print(np.mean(eta_MCMC))
                    print(eta_MCMC_high)
                    print(eta_MCMC_low)
                    axarr[0,0].errorbar(Tsat_sim[0],np.mean(eta_MCMC),yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc=color_list[model],markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)
                    if args.devPlot:
                        axarr[1,0].errorbar(Tsat_sim[0],per_error_MCMC_avg,yerr=np.array([[per_error_MCMC_high,per_error_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc=color_list[model],markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)
                except:
                    print('Could not find validation results for '+compound+' '+model+' irho = 0')

            axarr[0,0].text(label_position_list[compound][0],label_position_list[compound][1],label_list[compound])
            for iM, model in enumerate(model_list):
            
                root_path = compound+'/Gromacs/Saturation_Viscosity/'+model+tail_path_list[compound]+"_LINCS/"
                #print(root_path)

                #root_path = compound+'/Gromacs/Saturation_Viscosity/'+model
                #root_path = compound+'/Gromacs/'+model

                if os.path.isdir(root_path+'_N400/'): root_path += '_N400/'
                #elif os.path.isdir(root_path+'_N800/'): root_path += '_N800/'
                #elif os.path.isdir(root_path+'_N200/'): root_path += '_N200/'
                #elif os.path.isdir(root_path+'_N100/'): root_path += '_N100/'
                elif os.path.isdir(root_path): root_path += '/'
                else: 
                    print('No file path for '+compound+' with '+model)
                 
                if iM == 0: # os.path.exists(root_path+'SaturatedSettings.txt'): # and not ('sim_sat' in locals()): # iM == 0: #Only load the simulation temperatures and the RP/TDE data once
                
                    try:

                        sim_sat = np.loadtxt(root_path+'SaturatedSettings.txt',skiprows=1)

                    except:

                        sim_sat = np.loadtxt(root_path+'SaturationSettings.txt',skiprows=1)
                    
                    try:
                        Tsat_sim = sim_sat[:,2]
                    except:
                        Tsat_sim = np.array([sim_sat[2]])
                        
                    try:
                        #print("Searched for REFPROP at "+root_path+"REFPROP_eta_sat.txt")
                        RP_Tsat_eta = np.loadtxt(root_path+'REFPROP_eta_sat.txt',skiprows=1)
                        #print(1)
                        RP_Tsat = RP_Tsat_eta[:,0]
                        #print(2)
                        RP_eta_sat = RP_Tsat_eta[:,1] / 1000. #[Pa-s]
                        #print(3)
                        RP_eta_sim = np.interp(Tsat_sim,RP_Tsat,RP_eta_sat)
                        #print(4)

                        #print(5)
                        if args.uncer and REF_uncer[compound] >= .05:
                            #print(6)
                            axarr[0,0].plot(RP_Tsat,np.log10(RP_eta_sat*(1+REF_uncer[compound])) + shifts[compound],'k:',linewidth=1)
                            #print(7)
                            axarr[0,0].plot(RP_Tsat,np.log10(RP_eta_sat*(1-REF_uncer[compound])) + shifts[compound],'k:',linewidth=1)
                        RP_eta_sat = np.log10(RP_eta_sat) + shifts[compound]
                        axarr[0,0].plot(RP_Tsat,RP_eta_sat,'k-',linewidth=3)

                    except:

                        print('Could not find REFPROP values for '+compound)

                    try:
                        
                        TDE_Tsat_eta = np.loadtxt(root_path+'TDE_eta_sat.txt',skiprows=1)
                        #print('Tried TDE in '+root_path+'TDE_eta_sat.txt')
                        #print(TDE_Tsat_eta)
                        #print(TDE_Tsat_eta.shape)

                        try:
                            TDE_Tsat = TDE_Tsat_eta[:,0]
                            TDE_eta_sat = TDE_Tsat_eta[:,1] / 1000. #[Pa-s]
                            #print("Try successful for "+compound)
                        except:
                            TDE_Tsat = TDE_Tsat_eta[0]
                            TDE_eta_sat = TDE_Tsat_eta[1] / 1000. #[Pa-s]
                            #print("Except successful for "+compound)
                        TDE_eta_sat = np.log10(TDE_eta_sat)
                        axarr[0,0].plot(TDE_Tsat,TDE_eta_sat + shifts[compound],'kx',markersize=10,markeredgewidth=2,mfc='None')
                        #print("START NEW:"+compound)
                        #print(TDE_Tsat)
                        #print("BREAK")
                        #print(TDE_eta_sat)
                        if args.TDEuncer:
                            axarr[0,0].errorbar(TDE_Tsat,TDE_eta_sat, yerr=TDE_Uncer_sat,mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)  # NOTE I just did this
                        
                    except:

                        print('Could not find TDE data for '+compound)
                    
                    nrho = len(Tsat_sim)
                    irho0 = int(nrhomax[compound] - nrho)
    
                for irho in range(nrho):
                    
                    try:
                        eta_MCMC =  np.loadtxt(root_path+'GK_eta_boots_rho'+str(irho+irho0)) / 1000. #[Pa-s]
                        per_error_MCMC = (eta_MCMC-RP_eta_sim[irho])/RP_eta_sim[irho]*100.
                        eta_MCMC = np.log10(eta_MCMC) + shifts[compound]
                        #print("Searching "+root_path+'GK_eta_boots_rho'+str(irho+irho0))
                        eta_MCMC_sorted = np.sort(eta_MCMC)
                        i95 = int(0.025*len(eta_MCMC_sorted))
                        eta_MCMC_95 = eta_MCMC_sorted[i95:-i95]
                        eta_MCMC_avg = np.mean(eta_MCMC_95) #Only take average without the outliers          
                        eta_MCMC_low = eta_MCMC_avg - eta_MCMC_95[0]
                        eta_MCMC_high = eta_MCMC_95[-1] - eta_MCMC_avg
                        per_error_MCMC_sorted = np.sort(per_error_MCMC)
                        i95 = int(0.025*len(per_error_MCMC_sorted))
                        per_error_MCMC_plot = per_error_MCMC_sorted[i95:-i95]
                        per_error_MCMC_avg = np.mean(per_error_MCMC_plot) #Only take average without the outliers          
                        per_error_MCMC_low = per_error_MCMC_avg - per_error_MCMC_plot[0]
                        per_error_MCMC_high = per_error_MCMC_plot[-1] - per_error_MCMC_avg
                        print(Tsat_sim[0])
                        print(np.mean(eta_MCMC))
                        print(eta_MCMC_high)
                        print(eta_MCMC_low)
     
#                        axarr[0].errorbar(Tsat_sim[irho],np.mean(eta_MCMC),yerr=RP_eta_sim[irho]/100.*np.array([[per_error_MCMC_high,per_error_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)    
                        # use np.log10(np.mean...
                        axarr[0,0].errorbar(Tsat_sim[irho],np.mean(eta_MCMC),yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
                        if args.devPlot:
                            axarr[1,0].errorbar(Tsat_sim[irho],per_error_MCMC_avg,yerr=np.array([[per_error_MCMC_high,per_error_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)
                    except:
                        print('Could not find simulation results for '+compound+' '+model+' irho = '+str(irho))
        axarr[0,0].legend()
        axarr[0,0].set_xlabel('Temperature (K)')
        # CHANGE THIS 
        if args.settings == 1 or args.settings == 2 or args.settings == 3:
            axarr[0,0].set_ylabel(r'log$_{10}(\eta_{\rm liq}^{\rm sat}$/Pa-s) + $\Delta_{\eta}$')
        else:
            axarr[0,0].set_ylabel(r'log$_{10}(\eta_{\rm liq}^{\rm sat}$/Pa-s)')

        # Special axis cuts and labels
        if args.settings == 1:
            axarr[0,0].set_xlim([150, 530])
            axarr[0,0].text(450,-2.70,r'$\Delta_\eta$=1')
            axarr[0,0].text(450,-3.77,r'$\Delta_\eta$=0.5')
            axarr[0,0].text(395,-4.2,r'$\Delta_\eta$=0')
            axarr[0,0].text(430,-4.75,r'$\Delta_\eta$=-0.5')
        elif args.settings == 2:
            axarr[0,0].text(490,-2.70,r'$\Delta_\eta$=1.5')
            axarr[0,0].text(490,-4.0,r'$\Delta_\eta$=0')
            axarr[0,0].text(490,-5.05,r'$\Delta_\eta$=-1.5')
            axarr[0,0].set_xlim([200, 570])
        
        #axarr[0].set_xlim([0.99*Tsat_sim.min(),1.01*Tsat_sim.max()])
        #axarr[0].set_ylim([0.9*RP_eta_sim.min(),1.1*RP_eta_sim.max()])
        if args.devPlot: 
            axarr[1,0].set_xlabel('Temperature (K)')
            axarr[1,0].set_ylabel(r'$\left(\eta_{\rm sim} - \eta_{\rm REFPROP}\right)/\eta_{\rm REFPROP} \times 100$%')
            axarr[1,0].set_xlim(axarr[0,0].get_xlim())
            #if args.uncer:
            #    if args.settings == 0:
            #        axarr[1,0].plot([0,370],[5,5],'k--')
            #        axarr[1,0].plot([0,370],[-5,-5],'k--')
            #        axarr[1,0].text(200,10,"REFPROP Uncertainty")
            #    if args.settings == 1:
            #        axarr[1,0].plot([0,690],[10,10],'k--')
            #        axarr[1,0].plot([0,690],[-10,-10],'k--')
            #        axarr[1,0].text(400,7,"REFPROP Uncertainty")
            # Plot validation points if they exist
        fig.tight_layout()
        
        fig.savefig('Viscosity_compare_force_fields/compare_force_fields.pdf')

        plt.show()  
    
if __name__ == '__main__':
    '''
    python compare_force_fields.py --comp str --mod str --devPlot --uncer --TDEuncer --Settings
  
    '''

    main()
