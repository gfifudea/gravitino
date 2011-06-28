#!/usr/bin/env python
'''
Modules to run and analyse SPheno. Current implementation icludes:
* Gneration of general LesHouches dictionary.
* Change RVKAPPAIN and RVSNVEVIN from a six component vector
* Gneration of general LesHouches.in file.
* Filter Decays Block.
* ranlog: for generate random numbers in ranges of several orders
  of magnitude
* Funtions to analyse neutrino solutions
'''
import pyslha

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

def writeLHAinFile(xdict,lhinfile='LesHouches.in',universal=True,RP=False):
    '''To write LesHouches.in in the right order.'''
    LesHouches2={}
    LesHouches2['AMODSEL']=xdict['MODSEL']
    LesHouches2['BSMINPUTS']=xdict['SMINPUTS']
    LesHouches2['CMINPAR']=xdict['MINPAR']
    if not universal: LesHouches2['DEXTPAR']=xdict['EXTPAR']
    if not RP:
        LesHouches2['ERVSNVEVIN']=xdict['RVSNVEVIN']
        LesHouches2['FRVKAPPAIN']=xdict['RVKAPPAIN']
        
    LesHouches2['GSPhenoInput']=xdict['SPHENOINPUT']
    if not RP:
        LesHouches2['HNEUTRINOBOUNDSIN']=xdict['NEUTRINOBOUNDSIN']
        
    pyslha.writeSLHAFile(lhinfile,LesHouches2,{})


def filterParticle(p1,pdaug,sign=False):
    '''Find decays channels with PDG PI:
         pdaug
    for Particle class:
        p1
    Returns filtered Particles class:
        p2
    '''
    p2 = pyslha.Particle(p1.pid, p1.totalwidth, p1.mass)
    if sign:
        p2.decays = [d for d in p1.decays if pdaug in d.ids]
    else:
        p2.decays = [d for d in p1.decays if pdaug in map(abs, d.ids)]
    return p2

if __name__ == '__main__':
    print '''Check pyspheno_random.py and
    pysphenofull.py for examples of use'''
