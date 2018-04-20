from __future__ import division
import numpy as np 
import argparse

def comb_rule(xii,xjj,comb):
    if comb == 'geometric':
        
        xij = np.sqrt(xii*xjj)
        
    elif comb == 'arithmetic':
        
        xij = (xii + xjj)/2.
                 
    f = open('comb_rule','w')
    f.write(str(xij))
    f.close()
    
def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--comb",type=str,choices=['geometric','arithmetic'],help="choose which type of combining rule to use")
    parser.add_argument("-x","--xparam",type=float,nargs='+',help="Provide two parameters to combine" )
    args = parser.parse_args()
    comb_rule(args.xparam[0],args.xparam[1],args.comb)

if __name__ == '__main__':
    '''
    computes the combining rule
    '''
    
    main()