#!/usr/bin/env python
from pyspheno import * #also imports pyslha
import commands

if __name__ == '__main__':
    #Set the SPheno command, e.g: ./SPheno, or:
    sphenocmd='SPheno_intel'
    #Build SPhheno input file:
    LesHouches=buildLHAinFile() 
    #Be sure that the neutrino fit is switched on:
    LesHouches['SPHENOINPUT'].entries[90]=0
    #m12min=200;m12max=800
    #m0min=200;m0max=1000
    m12min=150;m12max=150
    m0min=1000;m0max=1000
    step=50
    for m12 in range(m12min,m12max+step,step): #must be integers
        for m0 in range(m0min,m0max+step,step): #must be integers
            LesHouches['MINPAR'].entries[1]=m0
            LesHouches['MINPAR'].entries[2]=m12
            writeLHAinFile(LesHouches,'LesHouches.in',RP=True)
            lsout=commands.getoutput(sphenocmd)
            spc,decays=pyslha.readSLHAFile('SPheno.spc')
            print commands.getoutput('slhaplot --br=10% --format png SPheno.spc')

    
#
