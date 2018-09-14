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

#a0 3.1015 3.7699 4.8865
#a1 7.6409 9.6459 12.9978
#a2 2.2294 2.2294 2.2294
#a3 -14.4432 -17.1165 -21.5858

a_AUA4_CH_CH2 = np.array([3.1015, 7.6409, 2.2294, -14.4432]) #[kJ/mol]
a_AUA4_CH_CH2_15 = np.array([3.7699,9.6459,2.2294,-17.1165]) #[kJ/mol]
a_AUA4_CH_CH2_40 = np.array([4.8865,12.9978,2.2294,-21.5858]) #[kJ/mol]

U_AUA4_CH_CH2 = U_plot_form1(a_AUA4_CH_CH2)-np.min(U_plot_form1(a_AUA4_CH_CH2))
U_AUA4_CH_CH2_15 = U_plot_form1(a_AUA4_CH_CH2_15)-np.min(U_plot_form1(a_AUA4_CH_CH2))
U_AUA4_CH_CH2_40 = U_plot_form1(a_AUA4_CH_CH2_40)-np.min(U_plot_form1(a_AUA4_CH_CH2))

fig = plt.figure(figsize=[12,12])

plt.plot(phi_plot,U_AUA4_CH_CH2*R_g,'g-',label='AUA4, CH-CH2')
plt.plot(phi_plot,U_AUA4_CH_CH2_15*R_g,'b:',label='AUA4m, CH-CH2, 15')
plt.plot(phi_plot,U_AUA4_CH_CH2_40*R_g,'r-.',label='AUA4m, CH-CH2, 40')
plt.plot(phi_plot,U_TraPPE_CH2_CH*R_g,'k--',label='TraPPE, CH2-CH')
plt.xlabel('Torsion (degrees)')
plt.ylabel(r'$U^{\rm tors}$ (K)')
plt.legend()
plt.show()

#a0 1.9176 2.1874 2.6639
#a1 5.7526 6.5618 7.9912
#a2 0.0 0.0 0.0
#a3 27.6703 28.7493 210.6552

a_AUA4_C_CH2 = np.array([1.9176,5.7526,0.0,-7.6703]) #[kJ/mol]
a_AUA4_C_CH2_15 = np.array([2.1874,6.5618,0.0,-8.7493]) #[kJ/mol]
a_AUA4_C_CH2_40 = np.array([2.6639,7.9912,0.0,-10.6552]) #[kJ/mol]

U_AUA4_C_CH2 = U_plot_form1(a_AUA4_C_CH2)
U_AUA4_C_CH2_15 = U_plot_form1(a_AUA4_C_CH2_15)
U_AUA4_C_CH2_40 = U_plot_form1(a_AUA4_C_CH2_40)

fig = plt.figure(figsize=[12,12])

plt.plot(phi_plot,U_AUA4_C_CH2,'g-',label='AUA4, C-CH2')
plt.plot(phi_plot,U_AUA4_C_CH2_15,'b:',label='AUA4m, C-CH2, 15')
plt.plot(phi_plot,U_AUA4_C_CH2_40,'r-.',label='AUA4m, C-CH2, 40')
plt.plot(phi_plot,U_TraPPE_CH2_C,'k--',label='TraPPE, CH2-C')
plt.xlabel('Torsion (degrees)')
plt.ylabel(r'$U^{\rm tors}$ (K)')
plt.legend()
plt.show()

c_AUA4_CH_CH2 = np.array([0.0,355.03,-68.19,791.35])

c0 = a_AUA4_CH_CH2[0] - a_AUA4_CH_CH2[1] + a_AUA4_CH_CH2[2] - a_AUA4_CH_CH2[3]  
c1 = a_AUA4_CH_CH2[1] + 3./4. * a_AUA4_CH_CH2[3] 
c2 = -a_AUA4_CH_CH2[2]/2.
c3 = a_AUA4_CH_CH2[3]/4.

print(c0/R_g)
print(c_AUA4_CH_CH2[0])

print(c1/R_g)
print(c_AUA4_CH_CH2[1])

print(c2/R_g)
print(c_AUA4_CH_CH2[2])

print(c3/R_g)
print(c_AUA4_CH_CH2[3])

#a_converted = 
#
#conversion_matrix = 