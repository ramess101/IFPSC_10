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

Mw_list = {'C3H8':44.096,'C4H10':58.122,'C8H18':114.23,'C12H26':170.33,'IC4H10':58.122,'IC8H18':114.23,'IC5H12':72.15,'3MPentane':86.2,'23DMButane':86.1754,'NEOC5H12':72.15}

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--comp",type=str,help="Specify the compound to analyze")
    parser.add_argument("--sat",action='store_true',help="If analyzing saturation conditions")
    parser.add_argument("--T293highP",action='store_true',help="If analyzing T=293 K and high pressures")
    parser.add_argument("--nrhomax",type=int,help="Specify the maximum number of rho values")
    parser.add_argument("-m","--mod",type=str,nargs='+',choices=['Potoff','TraPPE','TAMie','AUA4','Exp6','Mie_14','Mie_16'],help="Specify the model(s) to analyze")
    parser.add_argument("--devPlot",action='store_true',help="If deviation plots are desired")
    parser.add_argument("--uncer",action='store_true',help="If REFPROP uncertainty should be plotted")
    parser.add_argument("--TDEuncer",action='store_true',help="If TDE uncertainty should be plotted")
    args = parser.parse_args()

    try:
        comp = args.comp
        model_list = args.mod
    except:
        print('Please specify the compound and/or the models')
        return 0

    Mw = Mw_list[comp]
    
    if args.nrhomax:
        nrhomax = args.nrhomax
    else:
        nrhomax = 5

    color_list = {'Potoff':'r','TraPPE':'b','TAMie':'g','AUA4':'m','Mie_14':'b','Mie_16':'r'}
    line_list = {'C3H8':'-','C4H10':'--','C8H18':'-.','IC8H18': ':'}
    # NOTE these are likely incorrect for C3H8 and C4H10
    ref_uncer = {'C3H8':.06,'C4H10':.005,'C8H18':.02,'C12H26':.02,'IC4H10':.03,'IC5H12':.04,'IC6H14':.02,'IC8H18':.05,'NEOC5H12':.10,'3MPentane':.02,'23DMButane':.05,'C16H34':.02}
    symbol_list = {'Potoff':'o','TraPPE':'s','TAMie':'^','AUA4':'v','Mie_14':'s','Mie_16':'o'}
    rec_edge_list = {'C3H8':1000,'C4H10':200,'IC4H10':200}
    tail_path_list = {'C3H8':'_N400_LINCS/','C4H10':'_N400_LINCS/','C8H18':'_N400_LINCS/','C12H26':'_N400_LINCS/','IC4H10':'_N400_LINCS/','IC8H18':'_N400_LINCS/','C16H34':'_N200_LINCS','3MPentane':'_N400_LINCS/','NEOC5H12':'_N400_LINCS/','23DMButane':'_N400_LINCS/','IC5H12':'_N400_LINCS/'}
    left_text = {'C3H8':[661.9, -15],'IC4H10':[600,0],'C4H10':[570,6]}
    right_text = {'C3H8':[671.6, -5],'IC4H10':[600,0],'C4H10':[0,-3]}
   # if args.sat:
   #     REF_uncer = {'C3H8':.005,'C4H10':.005,'C8H18':.005,'C12H26':.005,'IC4H10':.03,'IC5H12':.04,'IC6H14':.02,'IC8H18':.05,'NEOC5H12':.10,'3MPentane':.02,'23DMButane':.05,'C16H34':.02}
    #label_list = {'C3H8':r'$n$-C$_3$','C4H10':r'$n$-C$_4$','C8H18':r'$n$-C$_8$','IC8H18':r'$i$-C$_8$'}
    #label_position_list = {'C3H8':[330,0.03],'C4H10':[400,0.03],'C8H18':[530,0.03],'IC8H18':[330,0.03]}
    #nrhomax = {'C3H8':10,'C4H10':5,'C8H18':5}
    
    for iM, model in enumerate(model_list):

        root_path = comp+'/Gromacs/'
        #if model == 'TraPPE': nrhomax = 3
        #if model == 'TAMie': nrhomax = 5
        if model == 'Mie_14' or model == 'Mie_16': tail_path_list[comp] = '_Bayesian_Viscosity_200_MCMC/'    
        if args.sat and args.T293highP:
            
                print('Must specify either saturation or high pressure conditions')
                return 0
        
        elif args.sat:
            
            print('This code does not work for saturation yet.')
            return 0
            
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
            RP_eta_sat = RP_Tsat_eta[:,1] /1000. #[Pa-s]
            RP_eta_sim = np.interp(Tsat_sim,RP_Tsat,RP_eta_sat)
            print("Boo! I am "+comp)
            print(RP_eta_sat)
        
            TDE_Tsat_eta = np.loadtxt('TDE_eta_sat.txt',skiprows=1)
        
            try:
                TDE_Tsat = TDE_Tsat_eta[:,0]
                TDE_eta_sat = TDE_Tsat_eta[:,1] /1000. #[Pa-s]
            except:
                TDE_Tsat = TDE_Tsat_eta[0]
                TDE_eta_sat = np.log10(TDE_Tsat_eta[1] /1000.) #[Pa-s] 
        
            fig,axarr = plt.subplots(ncols=2,nrows=1,figsize=[20,10])
    
            #axarr[0].set_yscale("log")
        
        elif args.T293highP:

            root_path = root_path+'T293highP_Viscosity/'+model+tail_path_list[comp]
        
            ### At some point saturated and high pressure should have different names here
        
            sim_sat = np.loadtxt(root_path+'T293highPSettings.txt',skiprows=1)
        
            try:
                T_sim = sim_sat[:,2]
                L_sim = sim_sat[:,1]
                N_sim = sim_sat[:,0]
            except:
                T_sim = np.array([sim_sat[2]])
                L_sim = np.array([sim_sat[1]])
                N_sim = np.array([sim_sat[0]])
            T_sim = T_sim[:nrhomax]
            L_sim = L_sim[:nrhomax]
            N_sim=N_sim[:nrhomax]           

            nrho = len(T_sim)
            print("nrho: "+str(nrho))
            press_sim = np.loadtxt(root_path+'press_all_log') #[bar]
            press_sim = press_sim[:60*nrhomax] 
            press_sim /= 10. #[MPa]
            print("length: "+str(len(press_sim)))
            print("60*nrho: "+str(60*nrhomax))

            press_sim = press_sim.reshape([nrho,int(len(press_sim)/nrho)])
        
            upress_sim = 1.96*np.std(press_sim,axis=1)
            press_sim = np.mean(press_sim,axis=1)

            try:

                L_sim = np.loadtxt(root_path+'Lbox_all') #[bar]
                L_sim = L_sim[:60*nrhomax]
        
                L_sim = L_sim.reshape([nrho,int(len(L_sim)/nrho)])
        
                uL_sim = 1.96*np.std(L_sim,axis=1)
                L_sim = np.mean(L_sim,axis=1)
            except:
                print('No Lbox_all file')
        
            rho_sim = Lbox_to_rho(L_sim,N_sim,Mw)
        
            RP_T293highP_eta = np.loadtxt(root_path+'REFPROP_eta_T293highP.txt',skiprows=1)
            #print("Boo! I am "+comp)
            #print(RP_T293highP_eta)
            RP_press = RP_T293highP_eta[:,1]
            RP_rho = RP_T293highP_eta[:,2]
            RP_eta_T293highP = RP_T293highP_eta[:,3] / 1000.
            RP_eta_rho_sim = np.interp(rho_sim,RP_rho,RP_eta_T293highP)
            RP_eta_press_sim = np.interp(press_sim,RP_press,RP_eta_T293highP)

            RP_uncer_high = np.log10(RP_eta_T293highP*(1+ref_uncer[comp]))
            RP_uncer_low = np.log10(RP_eta_T293highP*(1-ref_uncer[comp]))

            RP_eta_rho_sim = np.log10(RP_eta_rho_sim)
            RP_eta_press_sim = np.log10(RP_eta_press_sim)
            RP_eta_T293highP = np.log10(RP_eta_T293highP)
            end1 = RP_eta_rho_sim[-1]  # Where to stop plotting REFPROP uncertainties
            end2 = RP_eta_press_sim[-1]
            # Special cut offs necessary for properly plotting refprop extrapolation
            if comp in rec_edge_list.keys():
                cutoff = rec_edge_list[comp]
                RPD_eta_low = RP_eta_T293highP[RP_press < cutoff]
                RPD_press_low = RP_press[RP_press < cutoff]
                RPD_rho_low = RP_rho[RP_press < cutoff]
                RPD_eta_high = RP_eta_T293highP[RP_press >= cutoff]
                RPD_press_high = RP_press[RP_press >= cutoff]
                RPD_rho_high = RP_rho[RP_press >= cutoff]
                RP_uncer_high = RP_uncer_high[RP_press < cutoff]
                RP_uncer_low = RP_uncer_low[RP_press < cutoff]
                end1 = RPD_rho_high[1]  # Revised where to stop plotting REFPROP uncertainties
                end2 = RPD_press_high[1]

            TDE_TPeta = np.loadtxt(root_path+'TDE_eta_TP.txt',skiprows=1)
            TDE_press = TDE_TPeta[:,1]/1000. #[MPa]
            TDE_eta = np.log10(TDE_TPeta[:,2])  # [Pa-s] *1000. #[cP]
            TDE_Temp = TDE_TPeta[:,0]
            # TDE_Uncer = TDE_TPeta[:,3]  # Uncertainty in Pa-s

            TDE_press_T293 = TDE_press[TDE_Temp>291]
            TDE_eta_T293 = TDE_eta[TDE_Temp>291]
            TDE_Temp_T293 = TDE_Temp[TDE_Temp>291]
            # TDE_Uncer_T293 = TDE_Uncer[TDE_Temp>291]
            TDE_press_T293 = TDE_press_T293[TDE_Temp_T293<299]
            TDE_eta_T293 = TDE_eta_T293[TDE_Temp_T293<299]
            # TDE_Uncer_T293 = TDE_Uncer_T293[TDE_Temp_T293<299]
            TDE_Temp_T293 = TDE_Temp_T293[TDE_Temp_T293<299]

            if iM == 0:
                
                if args.devPlot:

                    fig,axarr = plt.subplots(ncols=2,nrows=2,figsize=[20,20])
   
                    #axarr[0,0].set_yscale("log")
                    #axarr[0,1].set_yscale("log")

                else:

                    fig,axarr = plt.subplots(ncols=2,nrows=1,figsize=[20,10])
   
                    #axarr[0].set_yscale("log")
                    #axarr[1].set_yscale("log")
    
        irho0 = int(nrhomax - nrho)
        
        for irho in range(nrho):
        
            try:
                eta_MCMC =  np.loadtxt(root_path+'GK_eta_boots_rho'+str(irho+irho0)) / 1000. #[Pa-s]
                temp_RP_eta = np.power(10,RP_eta_rho_sim[irho])
                per_error_MCMC = (eta_MCMC-temp_RP_eta)/temp_RP_eta*100.
                print("RP_eta_rho_sim")
                print(RP_eta_rho_sim)
                eta_MCMC_med = np.median(eta_MCMC)
                eta_MCMC = eta_MCMC[eta_MCMC > 0.75 * eta_MCMC_med]
                eta_MCMC = eta_MCMC[eta_MCMC < 1.25 * eta_MCMC_med]
                eta_MCMC = np.log10(eta_MCMC)
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
                    
                    if rho_sim[irho] > 200.: #Check to make sure this is a liquid density
                        per_error_MCMC_sorted = np.sort(per_error_MCMC)
                        i95 = int(0.025*len(per_error_MCMC_sorted))
                        per_error_MCMC_plot = per_error_MCMC_sorted[i95:-i95]
                        per_error_MCMC_avg = np.mean(per_error_MCMC_plot) #Only take average without the outliers          
                        per_error_MCMC_low = per_error_MCMC_avg - per_error_MCMC_plot[0]
                        per_error_MCMC_high = per_error_MCMC_plot[-1] - per_error_MCMC_avg 

                        if args.devPlot:

                            axarr[0,0].errorbar(rho_sim[irho],np.mean(eta_MCMC),yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
                
                            if rho_sim[irho] < RP_rho.max():
                                axarr[1,0].errorbar(rho_sim[irho],per_error_MCMC_avg,yerr=np.array([[per_error_MCMC_high,per_error_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)

                        else:

                            axarr[0].errorbar(rho_sim[irho],np.mean(eta_MCMC),yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
                
                        RP_eta_press_temp = np.power(10,RP_eta_press_sim[irho])
                        temp_MCMC = np.power(10,eta_MCMC) 
                        per_error_MCMC = (temp_MCMC-RP_eta_press_temp)/RP_eta_press_temp*100.
                        per_error_MCMC_sorted = np.sort(per_error_MCMC)
                        i95 = int(0.025*len(per_error_MCMC_sorted))
                        per_error_MCMC_plot = per_error_MCMC_sorted[i95:-i95]
                        per_error_MCMC_avg = np.mean(per_error_MCMC_plot) #Only take average without the outliers          
                        per_error_MCMC_low = per_error_MCMC_avg - per_error_MCMC_plot[0]
                        per_error_MCMC_high = per_error_MCMC_plot[-1] - per_error_MCMC_avg 

                        if args.devPlot:

                            axarr[0,1].errorbar(press_sim[irho],np.mean(eta_MCMC),xerr=upress_sim[irho],yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
                
                            if press_sim[irho] < RP_press.max():
                                axarr[1,1].errorbar(press_sim[irho],per_error_MCMC_avg,xerr=upress_sim[irho],yerr=np.array([[per_error_MCMC_high,per_error_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)
        
                        else:

                            axarr[1].errorbar(press_sim[irho],np.mean(eta_MCMC),xerr=upress_sim[irho],yerr=np.array([[eta_MCMC_high,eta_MCMC_low]]).T,fmt=color_list[model]+symbol_list[model],mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)                                         
       
            except:
             
                pass
    
        if args.sat and iM == len(model_list):
    
            axarr[0].plot(TDE_Tsat,TDE_eta_sat,'kx',markersize=10,markeredgewidth=2,mfc='None',label='TDE Exp. Data')
            axarr[0].plot(RP_Tsat,RP_eta_sat,'k-',linewidth=3,label='REFPROP')
            axarr[0].legend()
        
            axarr[0].set_xlabel('Temperature (K)')
            axarr[0].set_ylabel(r'log$_{10}$($\eta_{\rm liq}^{\rm sat}$/Pa-s)')

        
            #axarr[0].set_xlim([0.99*Tsat_sim.min(),1.01*Tsat_sim.max()])
            #axarr[0].set_ylim([0.9*RP_eta_sim.min(),1.1*RP_eta_sim.max()])
         
            axarr[1].set_xlabel('Temperature (K)')
            axarr[1].set_ylabel(r'$\left(\eta^{\rm sat}_{\rm sim} - \eta^{\rm sat}_{\rm REFPROP}\right)/\eta^{\rm sat}_{\rm REFPROP} \times 100$%')
            axarr[1].legend()   
        
            fig.tight_layout()
    
            fig.savefig('Viscosity_compare_force_fields/compare_TDE_REFPROP_sat_'+comp+'.pdf')
    
            plt.show()
        
        elif args.T293highP and iM == len(model_list)-1 and args.devPlot:


            for model in model_list: axarr[0,1].plot([],[],color_list[model]+symbol_list[model],markersize=10,markeredgewidth=2,mfc='None',label=model)
         
            if comp in rec_edge_list.keys():
                axarr[0,0].plot(RPD_rho_low,RPD_eta_low,'k-',linewidth=3,label='REFPROP')
                axarr[0,0].plot(RPD_rho_high,RPD_eta_high,'k--',linewidth=3,label='REFPROP Extrapolation')
            else:
                axarr[0,0].plot(RP_rho,RP_eta_T293highP,'k-',linewidth=3,label='REFPROP')

            if args.uncer:
                axarr[1,0].plot([0,end1],np.full(2,ref_uncer[comp] * 100),'k--')
                axarr[1,0].plot([0,end1],np.full(2,-ref_uncer[comp] * 100),'k--')
                axarr[1,1].plot([0,end2],np.full(2,ref_uncer[comp] * 100),'k--')
                axarr[1,1].plot([0,end2],np.full(2,-ref_uncer[comp] * 100),'k--')
                axarr[1,0].text(left_text[comp][0],left_text[comp][1])
                axarr[1,1].text(right_text[comp][0],right_text[comp][1])

            axarr[0,0].set_xlabel(r'$\rho$ (kg/m$^3$)')
            axarr[0,0].set_ylabel(r'log$_{10}$($\eta_{\rm liq}^{\rm comp}$/Pa-s)')
         
            axarr[1,0].set_xlabel(r'$\rho$ (kg/m$^3$)')
            axarr[1,0].set_ylabel(r'$\left(\eta_{\rm sim} - \eta_{\rm REFPROP}\right)/\eta_{\rm REFPROP} \times 100$%')
            axarr[0,1].plot(TDE_press_T293,TDE_eta_T293,'kx',markersize=10,markeredgewidth=2,mfc='None',label='TDE Exp. Data')
            if args.TDEuncer:
                axarr[0,1].errorbar(TDE_press_T293,TDE_eta_T293, yerr=TDE_Uncer_T293,fmt='r',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)  # NOTE I just did this
            if comp in rec_edge_list.keys():
                axarr[0,1].plot(RPD_press_low,RPD_eta_low,'k-',linewidth=3,label='REFPROP')
                axarr[0,1].plot(RPD_press_high,RPD_eta_high,'k--',linewidth=3,label='REFPROP Extrapolation')
            else:
                axarr[0,1].plot(RP_press,RP_eta_T293highP,'k-',linewidth=3,label='REFPROP')
            axarr[0,1].legend()

            # Plot REFPROP uncertainty
            if args.uncer and ref_uncer[comp] >= .05:
                if comp in rec_edge_list.keys():
                    axarr[0,0].plot(RPD_rho_low,RP_uncer_high,'k--',linewidth=1)
                    axarr[0,0].plot(RPD_rho_low,RP_uncer_low,'k--',linewidth=1)
                    axarr[0,1].plot(RPD_press_low,RP_uncer_high,'k--',linewidth=1)
                    axarr[0,1].plot(RPD_press_low,RP_uncer_low,'k--',linewidth=1)
                else:
                    axarr[0,0].plot(RP_rho,RP_uncer_high,'k--')
                    axarr[0,0].plot(RP_rho,RP_uncer_low,'k--')
                    axarr[1,0].plot(RP_press,RP_uncer_low,'k--')
                    axarr[1,0].plot(RP_press,RP_uncer_high,'k--')

            axarr[0,1].set_xlabel(r'$P$ (MPa)')
            axarr[0,1].set_ylabel(r'log$_{10}$($\eta_{\rm liq}^{\rm comp}$/Pa-s)')

         
            axarr[1,1].set_xlabel(r'$P$ (MPa)')
            axarr[1,1].set_ylabel(r'$\left(\eta_{\rm sim} - \eta_{\rm REFPROP}\right)/\eta_{\rm REFPROP} \times 100$%')

            axarr[1,0].set_xlim(axarr[0,0].get_xlim())
            axarr[1,1].set_xlim(axarr[0,1].get_xlim())

            
        
            fig.tight_layout()
    
            fig.savefig('Viscosity_compare_force_fields/compare_REFPROP_T293highP_'+comp+'.pdf')
    
            plt.show()

        elif args.T293highP and iM == len(model_list)-1:

            for model in model_list: axarr[1].plot([],[],color_list[model]+symbol_list[model],markersize=10,markeredgewidth=2,mfc='None',label=model)
         
            if comp in rec_edge_list.keys():
                axarr[0].plot(RPD_rho_low,RPD_eta_low,'k-',linewidth=3,label='REFPROP')
                axarr[0].plot(RPD_rho_high,RPD_eta_high,'k--',linewidth=3,label='REFPROP Extrapolation')
                if args.uncer and ref_uncer[comp] >= .05:
                  axarr[0].plot(RPD_rho_low,RP_uncer_high,'k--')
                  axarr[0].plot(RPD_rho_low,RP_uncer_low,'k--')
            else:
                axarr[0].plot(RP_rho,RP_eta_T293highP,'k-',linewidth=3,label='REFPROP')
                if args.uncer and ref_uncer[comp] >= .05:
                    axarr[0].plot(RP_rho,RP_uncer_high,'k--',linewidth=1)
                    axarr[0].plot(RP_rho,RP_uncer_low,'k--',linewidth=1)
            #axarr[0].legend()
            if args.TDEuncer:
                axarr[1].errorbar(TDE_press_T293,TDE_eta_T293, yerr=TDE_Uncer_T293,fmt='r',mfc='None',markersize=10,capsize=5,solid_capstyle='butt',markeredgewidth=2)  # NOTE I just did this
        
            axarr[0].set_xlabel(r'$\rho$ (kg/m$^3$)')
            axarr[0].set_ylabel(r'log$_{10}$($\eta_{\rm liq}^{\rm comp}$/Pa-s)')
         
            axarr[1].plot(TDE_press_T293,TDE_eta_T293,'kx',markersize=10,markeredgewidth=2,mfc='None',label='TDE Exp. Data')
            if comp in rec_edge_list.keys():
                axarr[1].plot(RPD_press_low,RPD_eta_low,'k-',linewidth=3,label='REFPROP')
                axarr[1].plot(RPD_press_high,RPD_eta_high,'k--',linewidth=3,label='REFPROP Extrapolation')
                if args.uncer and ref_uncer[comp] >= .05:
                  axarr[1].plot(RPD_press_low,RP_uncer_high,'k--')
                  axarr[1].plot(RPD_press_low,RP_uncer_low,'k--')
            else:
                axarr[1].plot(RP_press,RP_eta_T293highP,'k-',linewidth=3,label='REFPROP')
                if args.uncer and ref_uncer[comp] >= .05:
                    axarr[1].plot(RP_press,RP_uncer_high,'k--',linewidth=1)
                    axarr[1].plot(RP_press,RP_uncer_low,'k--',linewidth=1)

                
            axarr[1].legend()

            axarr[1].set_xlabel(r'$P$ (MPa)')
            axarr[1].set_ylabel(r'log$_{10}$($\eta_{\rm liq}^{\rm comp}$/Pa-s)')
                 
            fig.tight_layout()
    
            fig.savefig('Viscosity_compare_force_fields/compare_REFPROP_T293highP_'+comp+'.pdf')
    
            plt.show()

if __name__ == '__main__':
    '''
    python compare_TDE_REFPROP_force_fields.py --comp str --nrhomax int --sat --T293highP --mod str --devPlot --uncer --TDEuncer
  
    '''

    main()   
