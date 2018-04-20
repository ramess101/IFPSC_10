"""
Creates tabulated .xvg files for gromacs to use
"""
from __future__ import division
import numpy as np 
import argparse

kb = 1.38064852e-26 #[kJ/K]
NA = 6.02214e23 #[1/mol]

def create_tab(nrep,natt=6.,ncoul=1.):
    r = np.arange(0,2.401,0.0005)
    f = open('tab_it.xvg','w')
    rcutin = 0.05

    for ri in r:
        if ri < rcutin:
            U1 = rcutin**-ncoul
            F1 = ncoul * rcutin ** -(ncoul + 1.)
            U2 = -rcutin**-natt
            F2 = -natt * rcutin ** -(natt + 1.)
            U3 = rcutin**-nrep
            F3 = nrep * rcutin ** -(nrep + 1.)
        else:
            U1 = ri**-ncoul
            F1 = ncoul * ri ** -(ncoul + 1.)
            U2 = -ri**-natt
            F2 = -natt * ri ** -(natt + 1.)
            U3 = ri**-nrep
            F3 = nrep * ri ** -(nrep + 1.)
            
        f.write(str(ri)+'\t'+str(U1)+'\t'+str(F1)+'\t'+str(U2)+'\t'+str(F2)+'\t'+str(U3)+'\t'+str(F3)+'\n')
    f.close()
    
def convert_eps_sig_C6_Clam(eps,sig,lam,n=6.,print_Cit=True):
    Ncoef = lam/(lam-n)*(lam/n)**(n/(lam-n))
    eps *= kb * NA
    C6 = Ncoef * eps * sig ** n
    Clam = Ncoef * eps * sig ** lam
    
    if print_Cit:
    
        f = open('C6_it','w')
        f.write(str(C6))
        f.close()
        
        f = open('Clam_it','w')
        f.write(str(Clam))
        f.close()
        
    else:
        
        return C6, Clam

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-l","--lam",type=float,help="Set the value for lambda")
    parser.add_argument("-e","--epsilon",type=float,help="Set the value for epsilon (K)")
    parser.add_argument("-s","--sigma",type=float,help="Set the value for sigma (nm)")
    parser.add_argument("-LJ","--LennardJones",help="Flag if using LJ model",action="store_true")
    args = parser.parse_args()
    if args.lam:
        create_tab(args.lam)
        if args.epsilon and args.sigma:
            convert_eps_sig_C6_Clam(args.epsilon,args.sigma,args.lam)
    else:
        if args.LennardJones:
            if args.epsilon and args.sigma:
                convert_eps_sig_C6_Clam(args.epsilon,args.sigma,12.)
        else:
            print('Please specify a value for lambda if using Mie potential. Or specify model type as LJ.')

if __name__ == '__main__':
    
    main()
