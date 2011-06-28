#!/usr/bin/env python
'''
* ranlog: for generate random numbers in ranges of several orders
  of magnitude
* Funtions to analyse neutrino solutions
'''
from pyspheno import * #includes pyslha
import numpy as np
import scipy.optimize 
import commands
import sys

def changeLHAinFile(x,xdict):
    """Change specfic entries of the xdict dictionary:
    """
    for i in range(1,4):
        xdict['RVKAPPAIN'].entries[i]=x[i-1] # epsilon_i
        xdict['RVSNVEVIN'].entries[i]=x[i+2] # v_L_i

def oscilation(spcfile):
    """oscilation parameters"""
    slha=pyslha.readSLHAFile(spcfile)
    #neutrino parameters
    Delta2m32=slha[0]['SPHENORP'].entries[7]
    Delta2m21=slha[0]['SPHENORP'].entries[8]
    s223=np.sin(np.arctan(np.sqrt(slha[0]['SPHENORP'].entries[9])))**2
    s212=np.sin(np.arctan(np.sqrt(slha[0]['SPHENORP'].entries[10])))**2
    U13=slha[0]['SPHENORP'].entries[11]
    return Delta2m32,Delta2m21,s223,s212,U13


def check_slha(spcfile):
    '''Print neutrino parameters for spcfile
    input: A array (n,n)
           the output of matrixUm.
    '''
    dm21r=[7.05E-5,8.34E-5];dm23r=[2.07E-3,2.75E-3]
    s223r=[0.36,0.67];s212r=[0.25,0.37];U213=0.056
    (Delta2m32,Delta2m21,s223,s212,U13)=oscilation(spcfile)
    print ' %.2E < Delta m_{32}^2_{exp} < %.2E' %(dm23r[0],dm23r[1])
    print 'Delta m^2_{23}=%.2E' %(Delta2m32)
    print ' %.2E < Delta m_{21}^2_{exp} < %.2E' %(dm21r[0],dm21r[1])
    print 'Delta m^2_{21}=%.2E' %(Delta2m21)
    print 'tan(theta23)=U[1,2]/U[2,2]=>sin^2(theta23)=%.2f' %s223
    print ' %.2f < tan(theta23)_{exp} < %.2f' %(s223r[0],s223r[1])
    print 'tan(theta12)=U[0,1]/U[0,0]=>sin^2(theta12)=%.2f' %s212
    print ' %.2f < sin^2(theta12)_{exp} < %.2f' %(s212r[0],s212r[1])
    print 'U_13^2=%.3f' %(abs(U13**2))
    print 'U_13_{exp}<%.3f' %U213        

def chisq(x,xdict,sphenocmd):
    '''Function to be optimized
    input: x -> array
             x[0:3] 
             x[3:6]
             '''
    eps=x[0:3]
    vi=x[3:6]
    #change xdict dictionary:
    changeLHAinFile(x,xdict)
    #generate LesHouches file:
    writeLHAinFile(xdict,'LesHouches.in')
    #Run Spheno
    lsout=commands.getoutput(sphenocmd)
    if lsout:
        print lsout
        sys.exit()

    (Delta2m32,Delta2m21,s223,s212,U13)=oscilation('SPheno.spc')
    #constants going into the chisq:
    s212r=0.304;s223r=0.5
    U213=0.001
    m22=2.40;m21=7.65
    if U13**2<0.056:
        return ((s212-s212r)/0.07)**2+((s223-s223r)/0.165)**2\
           +((1E3*Delta2m32-m22)/0.4)**2+((1E5*Delta2m21-m21)/0.6)**2
    else:
        return 10.

def ranlog(vmin,vmax,dim=()):
    '''Random number generation uniformelly for
    a range expanding several orders of magnitude
    and random signs.'''
    return np.exp((np.log(vmax)-np.log(vmin))*np.random.uniform(0,1,dim)+np.log(vmin))*(-1)**np.random.random_integers(0,1,dim)

def searchmin(x0,xdict):
    '''Find the minimum'''
    #    return scipy.optimize.fmin_powell(chisq,x0,\
    #                xtol=1E-14,ftol=1E-14,full_output=1)[1]
    return scipy.optimize.fmin_powell(chisq,x0,args=(xdict,sphenocmd),\
                full_output=1)



def optloop(ifin,xdict,minimum=False,mu=100,vd=100):
    '''main Loop'''
    #np.random.seed(1)
    if minimum:
#        X0=np.random.uniform(-0.12,0.12,(ifin,6))
        eps=ranlog(1E-5,1,(ifin,3))
        Lam1=ranlog(1E-5,6E-02,(ifin,1))
        Lam=np.hstack([Lam1,ranlog(3E-2,7E-01,(ifin,2))])
        vL=(Lam-vd*eps)/mu
        X0=np.hstack([eps,vL])
    else:
        eps=ranlog(1E-5,1,(ifin,3))
        Lam1=ranlog(1E-5,6E-02,(ifin,1))
        Lam=np.hstack([Lam1,ranlog(3E-2,7E-01,(ifin,2))])
        vL=(Lam-vd*eps)/mu
        X0=np.hstack([eps,vL])

    if minimum:
        fmin=1E16
        for x0 in X0:
            sm=searchmin(x0,xdict)
            if sm[1] < fmin:
                fmin=sm[1]
                xmin=sm[0]
    else:
        fun=np.apply_along_axis(chisq,1,X0,xdict,sphenocmd)
        xmin=np.concatenate((X0,np.vstack(fun)),axis=1)

    return xmin
    
if __name__ == '__main__':
    m0=1000.
    m12=500.
    sphenocmd='SPheno_intel'
    LesHouches=buildLHAinFile() #see below
    LesHouches['MINPAR'].entries[1]=m0
    LesHouches['MINPAR'].entries[2]=m12

    minimum=False
    if minimum:
        ifin=1
        mu=100;vd=100 #not used at all
    else:
        ifin=1000

    #obtain vd and mu for msugra fixed parameters
    LesHouches['SPHENOINPUT'].entries[91]=1
    writeLHAinFile(LesHouches,'LesHouches.in')
    lsout=commands.getoutput(sphenocmd)
    if lsout:
        print lsout
        sys.exit()
        
    spc,decays=pyslha.readSLHAFile('SPheno.spc')
    vd=spc['SPHENORP'].entries[15]
    #mu=((Lamda-vd*epsilon)/vi)
    mu=(spc['SPHENORP'].entries.values()[1]\
        -vd*spc['RVKAPPA'].entries.values()[1])\
        /spc['RVSNVEV'].entries.values()[1]

    #Turn off neutrino fit
    LesHouches['SPHENOINPUT'].entries[91]=0

    B=optloop(ifin,LesHouches,minimum,mu,vd)
    x0=B[:,0:-1][B[:,-1].argmin()]
    #check neutrino fit Function
    print chisq(x0,LesHouches,sphenocmd)
    #print 
    check_slha('SPheno.spc')
    
#Fron RpParaneters.dat in SPhenoRP_intel/test
#3E-02 0.6 ! epsilon_1
#1.E-05 1. ! epsilon_2
#1.E-05 1. ! epsilon_3
#1.E-05 6.E-02 ! Lambda_1 = v_d epsilon_1 + mu v_L1
#3.E-02 7E-01 ! Lambda_2 = v_d epsilon_2 + mu v_L2
#3.E-02 7E-01 ! Lambda_3 = v_d epsilon_3 + mu v_L3


#
