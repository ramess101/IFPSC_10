# -*- coding: utf-8 -*-
"""
Class for quantifying uncertainty in torsional parameters

"""

from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import skewnorm

kb=1.38064852*10.**(-23.) #[J/K]
R_g = 8.3144598 / 1000. #[kJ/mol/K]

font = {'size' : '28'}
plt.rc('font',**font)

def U_tors_form1(phi,a_tors):
    phi_rad = phi * 2. * np.pi / 360. #[radians]
    phi_rad += np.pi # Must shift for IUPAC notation
    
    U_total = 0.
    
    for i, ai in enumerate(a_tors):
        
        U_total += ai*np.cos(phi_rad)**i
                            
    U_total /= R_g #[K]
                            
    return U_total
    
a_AUA4 = np.array([8.32556,17.70555,-2.51974,-30.03364,18.51362,16.34541,-37.32590,-14.43552,23.42457])
a_AUA4m_CH3 = np.array([12.46720,30.13030,-2.51974,-46.60000,18.51362,16.34541,-37.32590,-14.43552,23.42457])
a_AUA4m_CH2 = np.array([9.88053,22.37050,-2.51974,-36.25350,18.51362,16.34541,-37.32590,-14.43552,23.42457])

phi_plot = np.linspace(0,360,361)

U_plot_form1 = lambda a_tors: U_tors_form1(phi_plot,a_tors)

U_AUA4 = U_plot_form1(a_AUA4)
U_AUA4m_CH3 = U_plot_form1(a_AUA4m_CH3)
U_AUA4m_CH2 = U_plot_form1(a_AUA4m_CH2)

def U_tors_form2(phi,a_tors):
    phi_rad = phi * 2. * np.pi / 360. #[radians]
    
    U_total = 0.
    
    for i, ai in enumerate(a_tors):
        
        U_total += ai*(1.+(-1.)**(i+1.)*np.cos(i*phi_rad))
                            
    return U_total #[K]
    
a_TraPPE_CH2_CH2 = np.array([0.0,355.03,-68.19,791.32])
a_TraPPE_CH2_CH = np.array([-251.06,428.73,-111.85,441.27])
a_TraPPE_CH2_C = np.array([0.0,0.0,0.0,461.29])
a_TraPPE_CH_CH = a_TraPPE_CH2_CH.copy()

U_plot_form2 = lambda a_tors: U_tors_form2(phi_plot,a_tors)

U_TraPPE_CH2_CH2 = U_plot_form2(a_TraPPE_CH2_CH2)
U_TraPPE_CH2_CH = U_plot_form2(a_TraPPE_CH2_CH)
U_TraPPE_CH2_C = U_plot_form2(a_TraPPE_CH2_C)

def U_s(phi,A_s):
    
    phi_rad = phi * 2. * np.pi / 360. #[radians]
    phi_rad += np.pi # Must shift for IUPAC notation
    
    U_s = A_s * np.sin(3./2.*phi_rad)**2
    
    return U_s #[K]

U_plot_s = lambda A_s: U_s(phi_plot,A_s)

A_s_CH3 = np.max(U_AUA4m_CH3- U_AUA4)
A_s_CH2 = np.max(U_AUA4m_CH2- U_AUA4)

U_s_CH3 = U_plot_s(A_s_CH3)
U_s_CH2 = U_plot_s(A_s_CH2)

a_TraPPE_CH2_CH2_shift = a_TraPPE_CH2_CH2 + np.array([-A_s_CH2,0.,0.,A_s_CH2/2.])
U_TraPPE_CH2_CH2_shift = U_plot_form2(a_TraPPE_CH2_CH2_shift)

fig = plt.figure(figsize=[12,12])

plt.plot(phi_plot,U_AUA4,'k-',label='AUA4')
plt.plot(phi_plot,U_AUA4m_CH3,'r-',label='AUA4m_CH3')
plt.plot(phi_plot,U_AUA4m_CH2,'b-',label='AUA4m_CH2')
plt.plot(phi_plot,U_TraPPE_CH2_CH2,'k--',label='TraPPE_CH2_CH2')
plt.plot(phi_plot,U_TraPPE_CH2_CH,'r--',label='TraPPE_CH2_CH')
plt.plot(phi_plot,U_TraPPE_CH2_C,'b--',label='TraPPE_CH2_C')
plt.plot(phi_plot,U_s_CH3,'r:',label=r'U$_s$ CH3')
plt.plot(phi_plot,U_s_CH2,'b:',label=r'U$_s$ CH2')
plt.xlabel('Torsion (degrees)')
plt.ylabel(r'$U^{\rm tors}$ (K)')
plt.legend()
plt.show()

print(A_s_CH3,A_s_CH2)
print(r'Nieto-Draghi reported 40% increment in the terminal torsions, I compute:')
print(A_s_CH3/np.max(U_AUA4)*100.)
print(r'Nieto-Draghi reported 15% increment in the terminal torsions, I compute:')
print(A_s_CH2/np.max(U_AUA4)*100.)

U_s_shift = U_s_CH2

fig = plt.figure(figsize=[12,12])

plt.plot(phi_plot,U_TraPPE_CH2_CH2,'k-',label='TraPPE_CH2_CH2')
plt.plot(phi_plot,U_TraPPE_CH2_CH,'r-',label='TraPPE_CH2_CH')
plt.plot(phi_plot,U_TraPPE_CH2_C,'b-',label='TraPPE_CH2_C')
plt.plot(phi_plot,U_TraPPE_CH2_CH2+U_s_shift,'k--',label='TraPPE_CH2_CH2, Shifted')
plt.plot(phi_plot,U_TraPPE_CH2_CH+U_s_shift,'r--',label='TraPPE_CH2_CH, Shifted')
plt.plot(phi_plot,U_TraPPE_CH2_C+U_s_shift,'b--',label='TraPPE_CH2_C, Shifted')
plt.plot(phi_plot,U_s_shift,'m:',label=r'Shift')
plt.plot(phi_plot,U_TraPPE_CH2_CH2_shift,'y:',label='TraPPE_CH2_CH2, Shifted-alt')
plt.xlabel('Torsion (degrees)')
plt.ylabel(r'$U^{\rm tors}$ (K)')
plt.legend()
plt.show()

STD_CH2_CH2 = -0.15 * np.max(U_TraPPE_CH2_CH2) / 1.96
STD_CH2_CH = -0.15 * np.max(U_TraPPE_CH2_CH) / 1.96
STD_CH2_C = -0.15 * np.max(U_TraPPE_CH2_C) / 1.96                           
                         
N_MCMC = 20                        
                         
### Original way of generating normal data
#MCMC_CH2_CH2 = np.random.normal(0,STD_CH2_CH2,N_MCMC)
#MCMC_CH2_CH = np.random.normal(0,STD_CH2_CH,N_MCMC)
#MCMC_CH2_C = np.random.normal(0,STD_CH2_C,N_MCMC)   

### To generate a file with just shifted values, not sampled values
MCMC_CH2_CH2 = STD_CH2_CH2 * 1.96 * np.ones(N_MCMC)
MCMC_CH2_CH = STD_CH2_CH * 1.96 * np.ones(N_MCMC)
MCMC_CH2_C = STD_CH2_C * 1.96 * np.ones(N_MCMC) 
MCMC_CH_CH = MCMC_CH2_CH.copy()                    

### To generate normal data
#skew = .0
#loc = 0.
#### 15%
#upperU = 0.15
#### 40%
#upperU = 0.4

### To generate skewed data, with -15% and +80% as 95% confidence:
#skew = 1000.
#loc=-0.38
#upperU = 0.8

### To generate skewed data, with -40% and +80% as 95% confidence:
#skew=1.
#loc=0.
#upperU = 0.8

## To generate skewed data, with -15% and +40% as 95% confidence:
skew = 1.25
loc=0.
upperU = 0.4
#

try:
    upper95=skewnorm.interval(0.95,skew,loc=loc)[1]
    lower95=skewnorm.interval(0.95,skew,loc=loc)[0]
    
    lowerU = upperU*lower95/upper95
    
    print('Upper percent error = '+str(upperU*100)+r'%, lower percent error = '+str(lowerU*100)+r'%.')
        
    STD_CH2_CH2 = upperU * np.max(U_TraPPE_CH2_CH2) / upper95
    STD_CH2_CH = upperU * np.max(U_TraPPE_CH2_CH) / upper95
    STD_CH2_C = upperU * np.max(U_TraPPE_CH2_C) / upper95
    STD_CH_CH = STD_CH2_CH
                            
    MCMC_CH2_CH2 = skewnorm.rvs(skew,loc=loc*STD_CH2_CH2,scale=STD_CH2_CH2,size=N_MCMC)
    MCMC_CH2_CH = skewnorm.rvs(skew,loc=loc*STD_CH2_CH,scale=STD_CH2_CH,size=N_MCMC)
    MCMC_CH2_C = skewnorm.rvs(skew,loc=loc*STD_CH2_C,scale=STD_CH2_C,size=N_MCMC)   
    MCMC_CH_CH = skewnorm.rvs(skew,loc=loc*STD_CH_CH,scale=STD_CH_CH,size=N_MCMC)

    print('skew = '+str(skew)+', loc = '+str(loc))
    print(skewnorm.interval(0.95,skew,loc=loc))
    print(skewnorm.mean(skew,loc=loc))
    print(skewnorm.median(skew,loc=loc))
    print(skewnorm.interval(0.95,skew,loc=loc)[1]/skewnorm.interval(0.95,skew,loc=loc)[0])

except:
    
    print('skew, loc, and upperU not specified')
    pass

plt.hist(MCMC_CH2_CH2,bins=50,color='k',alpha=0.2)
plt.hist(MCMC_CH2_CH,bins=50,color='r',alpha=0.2)
plt.hist(MCMC_CH2_C,bins=50,color='b',alpha=0.2)
plt.show()

#skew_range = np.linspace(0,50,50)
#loc_range = np.linspace(-10,0,50)
#
#lower95 = np.zeros([len(skew_range),len(loc_range)])
#upper95 = np.zeros([len(skew_range),len(loc_range)])
#skewmean = np.zeros([len(skew_range),len(loc_range)])
#
#for iskew, skew in enumerate(skew_range):
#    for iloc, loc in enumerate(loc_range):
#
#        skewnorm_interval = skewnorm.interval(0.95,skew,loc=loc)
#        skewnorm_mean = skewnorm.mean(skew,loc=loc)
#        
#        lower95[iskew,iloc] = skewnorm_interval[0]
#        upper95[iskew,iloc] = skewnorm_interval[1]
#        skewmean[iskew,iloc] = skewnorm_mean
#        
##        if skewnorm_mean < 0.1 and skewnorm_mean > -0.1:
##            
##            print('skew = '+str(skew)+', loc = '+str(loc))
##            print(skewnorm.interval(0.95,skew,loc=loc))
##            print(skewnorm.mean(skew,loc=loc))
##            
##            if skewnorm_interval[0] < -1.5 and skewnorm_interval[1] > 4:
##                
##                print('skew = '+str(skew)+', loc = '+str(loc))
##                print(skewnorm.interval(0.95,skew,loc=loc))
##                print(skewnorm.mean(skew,loc=loc))
#
#fig = plt.figure(figsize=[6,6])
#CS = plt.contour(loc_range,skew_range,skewmean)
#cbar = fig.colorbar(CS)
#cbar.add_lines(CS)
#plt.title('Mean')
#plt.xlabel('loc')
#plt.ylabel('skew')
#plt.show()
#
#fig = plt.figure(figsize=[6,6])
#CS = plt.contour(loc_range,skew_range,lower95)
#cbar = fig.colorbar(CS)
#cbar.add_lines(CS)
#plt.title('Lower 95')
#plt.xlabel('loc')
#plt.ylabel('skew')
#plt.show()
#
#fig = plt.figure(figsize=[6,6])
#CS = plt.contour(loc_range,skew_range,upper95)
#cbar = fig.colorbar(CS)
#cbar.add_lines(CS)
#plt.title('Upper 95')
#plt.xlabel('loc')
#plt.ylabel('skew')
#plt.show()

fig = plt.figure(figsize=[8,8])

for MCMC_i in MCMC_CH2_CH2:
    
    a_MCMC = a_TraPPE_CH2_CH2 + np.array([-MCMC_i,0.,0.,MCMC_i/2.])
    U_MCMC = U_plot_form2(a_MCMC)
    
    plt.plot(phi_plot,U_MCMC,'r-',alpha=0.05)

plt.plot(phi_plot,U_TraPPE_CH2_CH2,'k-',label='TraPPE')
plt.plot([],[],'r-',label='MCMC')
plt.xlabel('Torsion (degrees)')
plt.ylabel(r'$U^{\rm tors}$ (K)')
plt.title(r'CH$_i$-CH$_2$-CH$_2$-CH$_j$')
plt.legend()
plt.show()

fig = plt.figure(figsize=[8,8])

for MCMC_i in MCMC_CH2_CH:
    
    a_MCMC = a_TraPPE_CH2_CH + np.array([-MCMC_i,0.,0.,MCMC_i/2.])
    U_MCMC = U_plot_form2(a_MCMC)
    
    plt.plot(phi_plot,U_MCMC,'r-',alpha=0.05)

plt.plot(phi_plot,U_TraPPE_CH2_CH,'k-',label='TraPPE')
plt.plot([],[],'r-',label='MCMC')
plt.xlabel('Torsion (degrees)')
plt.ylabel(r'$U^{\rm tors}$ (K)')
plt.title(r'CH$_i$-CH$_2$-CH-CH$_j$')
plt.legend()
plt.show()

fig = plt.figure(figsize=[8,8])

for MCMC_i in MCMC_CH2_C:
    
    a_MCMC = a_TraPPE_CH2_C + np.array([-MCMC_i,0.,0.,MCMC_i/2.])
    U_MCMC = U_plot_form2(a_MCMC)
    
    plt.plot(phi_plot,U_MCMC,'r-',alpha=0.05)

plt.plot(phi_plot,U_TraPPE_CH2_C,'k-',label='TraPPE')
plt.plot([],[],'r-',label='MCMC')
plt.xlabel('Torsion (degrees)')
plt.ylabel(r'$U^{\rm tors}$ (K)')
plt.title(r'CH$_i$-CH$_2$-C-CH$_j$')
plt.legend()
plt.show()

def build_GROMACS_tors(a_tors):
    
    F = np.zeros(6) #Using F and C to be consistent with notation in .top files
    C = np.zeros(6)
    
    for i, ai in enumerate(a_tors):
        
        F[i] = ai
    
    C[0] = F[0] + F[1] + 2.*F[2] + F[3]
    C[1] = -F[1] + 3.*F[3]
    C[2] = -2.*F[2] + 8*F[4]
    C[3] = -4.*F[3]
    C[4] = -8.*F[4]
    
    for i, ci in enumerate(C):
        
        if ci == -0:
            
            C[i] = 0. #Not sure if GROMACS would like negative zero entries
    
    C *= R_g # Convert to GROMACS units (kJ/mol)

    return C


print('GROMACS values for CHX  CH2  CH2  CHY:')
print('8.39736  16.78632   1.13393  -26.31760   0.00000   0.00000') 
print('Computed values:')
print(build_GROMACS_tors(a_TraPPE_CH2_CH2))   

print('GROMACS values for CHX  CH2  CH  CHY:')
print('3.28629   7.44211   1.85995  -14.67569   0.00000   0.00000') 
print('Computed values:')
print(build_GROMACS_tors(a_TraPPE_CH2_CH))  

print('GROMACS values for CHX  CH2  C  CHY:')
print('3.83538  11.50613   0.00000  -15.34151   0.00000   0.00000') 
print('Computed values:')
print(build_GROMACS_tors(a_TraPPE_CH2_C))  

#a_MCMC_CH2_CH2 = np.zeros([len(MCMC_CH2_CH2),len(a_TraPPE_CH2_CH2)])
#a_MCMC_CH2_CH = np.zeros([len(MCMC_CH2_CH),len(a_TraPPE_CH2_CH)])
#a_MCMC_CH2_C = np.zeros([len(MCMC_CH2_C),len(a_TraPPE_CH2_C)])
#
#for ijk, (MCMC_i, MCMC_j, MCMC_k) in enumerate(zip(MCMC_CH2_CH2,MCMC_CH2_CH,MCMC_CH2_C)):
#    
#    a_MCMC_CH2_CH2[ijk,:] = a_TraPPE_CH2_CH2 + np.array([-MCMC_i,0.,0.,MCMC_i/2.])
#    a_MCMC_CH2_CH[ijk,:] = a_TraPPE_CH2_CH2 + np.array([-MCMC_j,0.,0.,MCMC_j/2.])
#    a_MCMC_CH2_C[ijk,:] = a_TraPPE_CH2_CH2 + np.array([-MCMC_k,0.,0.,MCMC_k/2.])
#
#GROMACS_MCMC_CH2_CH2 = build_GROMACS_tors(a_MCMC_CH2_CH2)
#GROMACS_MCMC_CH2_CH = build_GROMACS_tors(a_MCMC_CH2_CH)
#GROMACS_MCMC_CH2_C = build_GROMACS_tors(a_MCMC_CH2_C)

GROMACS_MCMC_CH2_CH2 = np.zeros([len(MCMC_CH2_CH2),6])
GROMACS_MCMC_CH2_CH = np.zeros([len(MCMC_CH2_CH),6])
GROMACS_MCMC_CH2_C = np.zeros([len(MCMC_CH2_C),6])
GROMACS_MCMC_CH_CH = np.zeros([len(MCMC_CH_CH),6])

for ijkl, (MCMC_i, MCMC_j, MCMC_k, MCMC_l) in enumerate(zip(MCMC_CH2_CH2,MCMC_CH2_CH,MCMC_CH2_C,MCMC_CH_CH)):
    
    a_MCMC_CH2_CH2 = a_TraPPE_CH2_CH2 + np.array([-MCMC_i,0.,0.,MCMC_i/2.])
    a_MCMC_CH2_CH = a_TraPPE_CH2_CH + np.array([-MCMC_j,0.,0.,MCMC_j/2.])
    a_MCMC_CH2_C = a_TraPPE_CH2_C + np.array([-MCMC_k,0.,0.,MCMC_k/2.])
    a_MCMC_CH_CH = a_TraPPE_CH_CH + np.array([-MCMC_l,0.,0.,MCMC_l/2.])
    
    GROMACS_MCMC_CH2_CH2[ijkl,:] = build_GROMACS_tors(a_MCMC_CH2_CH2)
    GROMACS_MCMC_CH2_CH[ijkl,:] = build_GROMACS_tors(a_MCMC_CH2_CH)
    GROMACS_MCMC_CH2_C[ijkl,:] = build_GROMACS_tors(a_MCMC_CH2_C)
    GROMACS_MCMC_CH_CH[ijkl,:] = build_GROMACS_tors(a_MCMC_CH_CH)

f = open('MCMC_CH2_CH2_tors','w')
#f.write('C0'+'\t'+'C1'+'\t'+'C2'+'\t'+'C3'+'\t'+'C4'+'\t'+'C5'+'\n')

g = open('MCMC_CH2_CH_tors','w')
#g.write('C0'+'\t'+'C1'+'\t'+'C2'+'\t'+'C3'+'\t'+'C4'+'\t'+'C5'+'\n')

h = open('MCMC_CH2_C_tors','w')
#h.write('C0'+'\t'+'C1'+'\t'+'C2'+'\t'+'C3'+'\t'+'C4'+'\t'+'C5'+'\n')

i = open('MCMC_CH_CH_tors','w')
#i.write('C0'+'\t'+'C1'+'\t'+'C2'+'\t'+'C3'+'\t'+'C4'+'\t'+'C5'+'\n')

for MCMC_f,MCMC_g,MCMC_h,MCMC_i in zip(GROMACS_MCMC_CH2_CH2,GROMACS_MCMC_CH2_CH,GROMACS_MCMC_CH2_C,GROMACS_MCMC_CH_CH):
        
    for C_fj, C_gj, C_hj, C_ij in zip(MCMC_f,MCMC_g,MCMC_h,MCMC_i):
        
        f.write(str(C_fj)+'\t')
        g.write(str(C_gj)+'\t')
        h.write(str(C_hj)+'\t')
        i.write(str(C_ij)+'\t')
                
    f.write('\n')
    g.write('\n')
    h.write('\n')
    i.write('\n')    
    
f.close()
g.close()
h.close()
i.close()