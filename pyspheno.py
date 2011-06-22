#!/usr/bin/env python
import numpy as np
import scipy.optimize 
import commands
import pyslha
import sys

#DEBUG: get rid of global variable: LesHouches.
#Must have the MASS block to be read
#LesHouches,decaysin=pyslha.readSLHAFile('LesHouches_MASS.in')

def buildLHAinFile(universal=True):
    """Usage of pyslha.Block to generate one LesHouches.in file"""
    #Initialize dictionary of block clases
    xdict={}
    #define block class
    MODSEL=pyslha.Block('MODSEL') # Select model
    #Add entries to the class
    MODSEL.entries[1]=1 # mSUGRA
    #Add class to dictionary
    xdict['MODSEL']=MODSEL
    #====================
    SMINPUTS=pyslha.Block('SMINPUTS')   # Standard Model inputs
    SMINPUTS.entries[1]=1.279340E+02       # alpha_rm^-1(M_Z), MSbar, SM
    SMINPUTS.entries[2]=1.166390E-05       # G_F, Fermi constant
    SMINPUTS.entries[3]=1.172000E-01       # alpha_s(MZ) SM MSbar
    SMINPUTS.entries[4]=9.118760E+01       # Z-boson pole mass
    SMINPUTS.entries[5]=4.250000E+00       # m_b(mb) SM MSbar
    SMINPUTS.entries[6]=1.727000E+02       # m_top(pole)
    SMINPUTS.entries[7]=1.777000E+00       # m_tau(pole)
    xdict['SMINPUTS']=SMINPUTS
    #====================
    MINPAR=pyslha.Block('MINPAR') # Input parameters
    MINPAR.entries[1]=1000                # M0
    MINPAR.entries[2]=240                # m_1/2
    MINPAR.entries[3]=1.000000E+01       # tanb
    MINPAR.entries[4]=1                  # sign(mu)
    MINPAR.entries[5]=-1.000000E+02       # A0
    xdict['MINPAR']=MINPAR
    #====================
    if not universal:
        EXTPAR=pyslha.Block('EXTPAR')     # Gaugino masses
        EXTPAR.entries[1]=240      # M_1  240
        EXTPAR.entries[2]=240      # M_2  240 134 270
        EXTPAR.entries[3]=240      # M_3  240
        xdict['EXTPAR']=EXTPAR
    #====================
    RVSNVEVIN=pyslha.Block('RVSNVEVIN')   # sneutrino vevs at Q
    RVSNVEVIN.entries[1]=-6.12717708E-03  # v_L_1
    RVSNVEVIN.entries[2]= 6.60031840E-03  # v_L_2
    RVSNVEVIN.entries[3]=-6.34389921E-03  # v_L_3
    xdict['RVSNVEVIN']=RVSNVEVIN
    #====================
    RVKAPPAIN=pyslha.Block('RVKAPPAIN')  # bilinear RP parameters at Q
    RVKAPPAIN.entries[1]= 1.60784795E-01  # epsilon_1
    RVKAPPAIN.entries[2]=-1.67970270E-01  # epsilon_2
    RVKAPPAIN.entries[3]= 1.71230537E-01  # epsilon_3
    xdict['RVKAPPAIN']=RVKAPPAIN
    #====================
    SPHENOINPUT=pyslha.Block('SPHENOINPUT')  # SPheno specific input
    SPHENOINPUT.entries[ 1]= 0		    # error level
    SPHENOINPUT.entries[ 2]= 0		    # if =1, then SPA conventions are used
    SPHENOINPUT.entries[11]= 1		    # calculate branching ratios
    SPHENOINPUT.entries[12]= 1.00000000E-06 # write only branching ratios larger than this value
    SPHENOINPUT.entries[21]= 0		    # calculate cross section 
    SPHENOINPUT.entries[22]= 5.00000000E+02 # cms energy in GeV
    SPHENOINPUT.entries[23]= 0.00000000E+00 # polarisation of incoming e- beam
    SPHENOINPUT.entries[24]= 0.00000000E+00 # polarisation of incoming e+ beam
    SPHENOINPUT.entries[25]= 0		    # if 0 no ISR is calculated, if 1 ISR is caculated
    SPHENOINPUT.entries[26]= 1.00000000E-02 # write only cross sections larger than this value [fb]
    SPHENOINPUT.entries[31]=-1.00000000E+00 # m_GUT, if < 0 than it determined via g_1=g_2
    SPHENOINPUT.entries[32]= 0              # require strict unification g_1=g_2=g_3 if '1' is set
    SPHENOINPUT.entries[41]= 2.49520000E+00 # width of the Z-boson
    SPHENOINPUT.entries[42]= 2.11800000E+00 # width of the W-boson
    SPHENOINPUT.entries[90]= 1		    # if R-parity is added
    SPHENOINPUT.entries[91]= 1		    # Neutrino Fit on
    xdict['SPHENOINPUT']=SPHENOINPUT
    #====================
    NEUTRINOBOUNDSIN=pyslha.Block('NEUTRINOBOUNDSIN')  # The values are from 0808.2016 for 3sigma 
    NEUTRINOBOUNDSIN.entries[1 ]=2.07000000E-21	     # Delta m^2 atm, min
    NEUTRINOBOUNDSIN.entries[2 ]=2.75000000E-21	      # Delta m^2 atm, max
    NEUTRINOBOUNDSIN.entries[3 ]=5.62500000E-01	      # tan^2 atm, min
    NEUTRINOBOUNDSIN.entries[4 ]=2.03000000E+00	      # tan^2 atm, max
    NEUTRINOBOUNDSIN.entries[5 ]=7.05000000E-23	      # lower bound on solar
    NEUTRINOBOUNDSIN.entries[6 ]=8.34000000E-23	      # upper bound on solar
    NEUTRINOBOUNDSIN.entries[7 ]=3.33333300E-01	      # tan^2 sol, min
    NEUTRINOBOUNDSIN.entries[8 ]=5.87300000E-01	      # tan^2 sol, max
    NEUTRINOBOUNDSIN.entries[9 ]=0.00000000E+00	      # U^2e3,min
    NEUTRINOBOUNDSIN.entries[10]= 5.6000000E-02	      # U^2e3,max
    xdict['NEUTRINOBOUNDSIN']=NEUTRINOBOUNDSIN
    #====================
    return xdict 

def changeLHAinFile(x):
    """Change specfic entries of the global dictionary:
    LesHouches"""
    for i in range(1,4):
        LesHouches['RVKAPPAIN'].entries[i]=x[i-1] # epsilon_i
        LesHouches['RVSNVEVIN'].entries[i]=x[i+2] # v_L_i

def writeLHAinFile(xdict,lhinfile='LesHouches.in',universal=True):
    '''To write LesHouches.in in the right order.'''
    LesHouches2={}
    LesHouches2['AMODSEL']=xdict['MODSEL']
    LesHouches2['BSMINPUTS']=xdict['SMINPUTS']
    LesHouches2['CMINPAR']=xdict['MINPAR']
    if not universal: LesHouches2['DEXTPAR']=xdict['EXTPAR']
    LesHouches2['ERVSNVEVIN']=xdict['RVSNVEVIN']
    LesHouches2['FRVKAPPAIN']=xdict['RVKAPPAIN']
    LesHouches2['GSPhenoInput']=xdict['SPHENOINPUT']
    LesHouches2['HNEUTRINOBOUNDSIN']=xdict['NEUTRINOBOUNDSIN']
    pyslha.writeSLHAFile(lhinfile,LesHouches2,{})


def filterDecay(xdict,pid,pdaug):
    '''Build filtered pyslha decay dictionary from:
         xdict
       original pyslha dictionary, for particle number:
         pid
       decaying into a set of particles which includes
         pdaug
       between the daughters.
    '''
    prtclnda=pyslha.Particle(pid,xdict[pid].totalwidth,xdict[pid].mass)
    filterdict={}
    list1=[channel for channel in xdict[pid].decays if channel.nda==2 \
                and (abs(channel.ids[0])==pdaug or abs(channel.ids[1])==pdaug)]
    prtclnda.decays=list1+[channel for channel in xdict[pid].decays \
                        if channel.nda==3 and \
                      (abs(channel.ids[0])==pdaug or abs(channel.ids[1])==pdaug \
                         or abs(channel.ids[2])==pdaug)]
    filterdict[pid]=prtclnda
    return filterdict

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
    vi=x[3:6]
    #change LesHouches dictionary:
    changeLHAinFile(x)
    #generate LesHouches file:
    writeLHAinFile(LesHouches,'LesHouches.in')
    #Run Spheno
    lsout=commands.getoutput(sphenocmd)
    (Delta2m32,Delta2m21,s223,s212,U13)=oscilation('SPheno.spc')
    #constants going into the chisq:
    s212r=0.304;s223r=0.5
    U213=0.001
    m22=2.40E-3;m21=7.65E-5
    return np.abs(s212-s212r)/s212r+np.abs(s223-s223r)/s223r\
           +np.abs(Delta2m32-m22)/m22+np.abs(Delta2m21-m21)/m21\
           +np.abs(U13**2)#-U213)/U213

def ranlog(vmin,vmax,dim=()):
    '''Random number generation uniformelly for
    a range expanding several orders of magnitude'''
    return np.exp((np.log(vmax)-np.log(vmin))*np.random.uniform(0,1,dim)+np.log(vmin))*(-1)**np.random.random_integers(0,1,dim)

def searchmin(x0):
    '''Find the minimum'''
    #    return scipy.optimize.fmin_powell(chisq,x0,\
    #                xtol=1E-14,ftol=1E-14,full_output=1)[1]
    return scipy.optimize.fmin_powell(chisq,x0,\
                full_output=1)



def optloop(ifin=1,minimum=False,mu=100,vd=100):
    '''main Loop'''
    np.random.seed(1)
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
            sm=searchmin(x0)
            if sm[1] < fmin:
                fmin=sm[1]
                xmin=sm[0]
    else:
        argfmin=np.array([chisq(x0) for x0 in X0]).argmin()
        xmin=X0[argfmin]

    return xmin
    
if __name__ == '__main__':
    global sphenocmd,LesHouches
    sphenocmd='SPheno_intel'
    LesHouches=buildLHAinFile() #see below
    #to modify the dictionary of clases:
    LesHouches['SPHENOINPUT'].entries[91]=0
    minimum=True #False
    if minimum:
        ifin=1
        mu=100;vd=100 #not used at all
    else:
        ifin=1#100000

    #obtain vd and mu for msugra fixed parameters
    LesHouches['SPHENOINPUT'].entries[91]=1
    writeLHAinFile(LesHouches,'LesHouches.in')
    lsout=commands.getoutput(sphenocmd)
    spc,decays=pyslha.readSLHAFile('SPheno.spc')
    vd=spc['SPHENORP'].entries[15]
    #mu=((Lamda-vd*epsilon)/vi)
    mu=(spc['SPHENORP'].entries.values()[1]\
        -vd*spc['RVKAPPA'].entries.values()[1])\
        /spc['RVSNVEV'].entries.values()[1]

    #Turn off neutrino fit
    LesHouches['SPHENOINPUT'].entries[91]=0

    x0=optloop(ifin,minimum,mu,vd)        
    #check neutrino fit Function
    print chisq(x0)
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
