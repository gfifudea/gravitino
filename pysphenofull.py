#!/usr/bin/env python
import pyspheno
import commands
if __name__ == '__main__':
    global sphenocmd 
    sphenocmd='SPheno_intel'
    LesHouches=pyspheno.buildLHAinFile() #see below

    #Loop over initial x
    #obtain vd and mu for msugra fixed parameters
    LesHouches['SPHENOINPUT'].entries[91]=1
    #m12min=200;m12max=800
    #m0min=200;m0max=1000
    m12min=70;m12max=70
    m0min=100;m0max=100
    step=50
    for m12 in pyspheno.np.arange(m12min,m12max+step,step):
        for m0 in pyspheno.np.arange(m0min,m0max+step,step):
            LesHouches['MINPAR'].entries[1]=m0
            LesHouches['MINPAR'].entries[2]=m12
            pyspheno.writeLHAinFile(LesHouches,'LesHouches.in')
            lsout=commands.getoutput(sphenocmd)
            spc,decays=pyspheno.pyslha.readSLHAFile('SPheno.spc')
            chietau=pyspheno.filterDecay(pyspheno.filterDecay(decays,1000022,11),1000022,15)
            chietaubr=sum([x.br for x in chietau[1000022].decays])
            chimutau=pyspheno.filterDecay(pyspheno.filterDecay(decays,1000022,13),1000022,15)
            chimutaubr=sum([x.br for x in chimutau[1000022].decays])
            print 'solar mixing %d %d %g' %(m12,m0,chietaubr/chimutaubr)

    
#Fron RpParaneters.dat in SPhenoRP_intel/test
#3E-02 0.6 ! epsilon_1
#1.E-05 1. ! epsilon_2
#1.E-05 1. ! epsilon_3
#1.E-05 6.E-02 ! Lambda_1 = v_d epsilon_1 + mu v_L1
#3.E-02 7E-01 ! Lambda_2 = v_d epsilon_2 + mu v_L2
#3.E-02 7E-01 ! Lambda_3 = v_d epsilon_3 + mu v_L3


#
