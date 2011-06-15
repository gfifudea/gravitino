#!/usr/bin/env python
import numpy as np
import commands
import pyslha

global sphenocmd,LesHouches
sphenocmd='SPheno_intel'

def dictinputfile():
    """generate LesHouches.in"""
    #Initialize dictionary of block clases
    LesHouches={}
    #define block class
    MODSEL=pyslha.Block('MODSEL')
    #Add entrires to the class
    MODSEL.entries[1]=1
    #Add class to dictionary
    LesHouches['MODSEL']=MODSEL
    #====================
    SMINPUTS=pyslha.Block('SMINPUTS')   # Standard Model inputs
    SMINPUTS.entries[1]=   1.279340E+02       # alpha_rm^-1(M_Z), MSbar, SM
    SMINPUTS.entries[2]=   1.166390E-05       # G_F, Fermi constant
    SMINPUTS.entries[3]=   1.172000E-01       # alpha_s(MZ) SM MSbar
    SMINPUTS.entries[4]=   9.118760E+01       # Z-boson pole mass
    SMINPUTS.entries[5]=   4.250000E+00       # m_b(mb) SM MSbar
    SMINPUTS.entries[6]=   1.727000E+02       # m_top(pole)
    SMINPUTS.entries[7]=   1.777000E+00       # m_tau(pole)
    LesHouches['SMINPUTS']=SMINPUTS
    #====================
    return LesHouches 

#Set dictionary for input file as a global variable to avoid multiple calls
LesHouches=dictinputfile()

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

def chisq(x):
    '''Function to be optimized
    input: x -> array
             x[0:3] 
             x[3:6]
             '''
    eps=x[0:3]
    Lam=x[3:6]
    #generate LesHouches.in
    #.....
    #Run Spheno
    lsout=commands.getoutput(sphenocmd)
    (Delta2m32,Delta2m21,s223,s212,U13)=oscilation('SPheno.spc')
    #constants for chisq
    s212r=0.304;s223r=0.5
    U213=0.001
    m22=2.40E-3;m21=7.65E-5
    return np.abs(s212-s212r)/s212r+np.abs(s223-s223r)/s223r+np.abs(Delta2m32-m22)/m22+np.abs(Delta2m21-m21)/m21+np.abs(U13**2)#-U213)/U213


def sphenorndm(epsmin,epsmax,Lammin,Lammax,mu,vd,sign):
    """DEBUG"""
    eps=np.random.uniform(-1,1,3)
    Lambda=np.random.uniform(-1,1,3)
    vi=(Lambda-eps*vd)/mu

def random_search():
    """DEBUG"""
    epsmax=np.array([1,1,1])
    epsmin=np.array([-1,-1,-1])
    Lammax=np.array([1,1,1])
    Lammin=np.array([-1,-1,-1])
    sign=False
    sgn=lambda n: (-1)**np.random.random_integers(1,2,n)



if __name__ == '__main__':
    #Change some entry of LesHouches dictionary:
    LesHouches['SMINPUTS'].entries[1]=1.28E+02
    #Write SLHA file
    pyslha.writeSLHAFile('LesHouches_new.in',LesHouches,{})
    #obtain vd, and mu
    vd=200.
    mu=400.
    #Loop over initial x
    x=np.random.uniform(-1,1,3)
    #find the minumum
    #.....
    print chisq(x)



#3E-02 0.6 ! epsilon_1
#1.E-05 1. ! epsilon_2
#1.E-05 1. ! epsilon_3
#1.E-05 6.E-02 ! Lambda_1 = v_d epsilon_1 + mu v_L1
#3.E-02 7E-01 ! Lambda_2 = v_d epsilon_2 + mu v_L2
#3.E-02 7E-01 ! Lambda_3 = v_d epsilon_3 + mu v_L3


#
