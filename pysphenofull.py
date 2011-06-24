#!/usr/bin/env python
from pyspheno import * #also imports pyslha
import commands

if __name__ == '__main__':
    #Set the SPheno command, e.g: ./SPheno, or:
    sphenocmd='SPheno_intel'
    #Build SPhheno input file:
    LesHouches=buildLHAinFile() 
    #Be sure that the neutrino fit is switched on:
    LesHouches['SPHENOINPUT'].entries[91]=1
    #m12min=200;m12max=800
    #m0min=200;m0max=1000
    m12min=240;m12max=240
    m0min=1000;m0max=1000
    step=50
    for m12 in range(m12min,m12max+step,step): #must be integers
        for m0 in range(m0min,m0max+step,step): #must be integers
            LesHouches['MINPAR'].entries[1]=m0
            LesHouches['MINPAR'].entries[2]=m12
            writeLHAinFile(LesHouches,'LesHouches.in')
            lsout=commands.getoutput(sphenocmd)
            #decays is a dictionary containing the particle objects with keys
            #given by PDG particle identificaction numbers (PIN).
            #Check pyslha doc:
            spc,decays=pyslha.readSLHAFile('SPheno.spc')
            #Extract particle object with key PIN=1000022, neutralino, with
            #only decays channels having particle PIN=11 and PIN=15.
            chietau=filterParticle(filterParticle(decays[1000022],11),15)
            chietaubr=sum([x.br for x in chietau.decays])
            chimutau=filterParticle(filterParticle(decays[1000022],13),15)
            chimutaubr=sum([x.br for x in chimutau.decays])
            print 'solar mixing %d %d %g' %(m12,m0,chietaubr/chimutaubr)

    
#
