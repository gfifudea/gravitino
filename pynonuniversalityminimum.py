#!/usr/bin/env python
from pynonuniversality import *
import scipy.optimize 
def fmin(x,delta2,m0,m12):
    x=np.asarray(x)
    delta1=x[0]
    LesHouches=buildLHAinFile(universal=False) #see below
    LesHouches['MINPAR'].entries[1]=m0
    LesHouches['MINPAR'].entries[2]=m12
    LesHouches['EXTPAR'].entries[1]=m12*(1.+delta1)
    LesHouches['EXTPAR'].entries[2]=m12*(1.+delta2)
    LesHouches['EXTPAR'].entries[3]=m12
    LesHouches['SPHENOINPUT'].entries[91]=1
    writeLHAinFile(LesHouches,'LesHouches.in',universal=False)
    lsout=commands.getoutput(sphenocmd)
    if not lsout:
        spc,decays=pyslha.readSLHAFile(spcfile)
        #M=np.asarray([spc['EXTPAR'].entries[i] for i in range(4)])
        Ni1=np.asarray([spc['RVNMIX'].entries[i][1] for i in range(1,4)])
        Ni2=np.asarray([spc['RVNMIX'].entries[i][2] for i in range(1,4)])
        return np.dot(np.cos(thw)*Ni1+np.sin(thw)*Ni2,\
                              np.cos(thw)*Ni1+np.sin(thw)*Ni2)
    else:
        return 100

if __name__ == '__main__':
    '''Be sure to set sphenocmd:'''
    #sphenocmd='SPheno_intel'
    sphenocmd='SPheno3.1.4'
    spcfile='SPheno.spc'
    xx=[]
    m0=1000.
    m12=500.
    #Check universal point and get SM parameters:
    LesHouches=buildLHAinFile()
    LesHouches['MINPAR'].entries[1]=m0
    LesHouches['MINPAR'].entries[2]=m12
    LesHouches['SPHENOINPUT'].entries[91]=1
    writeLHAinFile(LesHouches,'LesHouches.in')
    lsout=commands.getoutput(sphenocmd)
    if not lsout:
        s2w,g2=smparameters(spcfile)
        thw=np.arcsin(np.sqrt(s2w))
    else:
        print "Wrong SUGRA point"
        sys.exit()


    delta1=1.
    delta2=0.
    U2sfnuuni=fmin([0],0.,m0,m12)
    x=np.array([delta1])
    print m0,m12,U2sfnuuni
    a=scipy.optimize.fmin_powell(fmin,x,args=(delta2,m0,m12),full_output=1)
    delta1=a[0]
    U2sfnumin=a[1]
    xx.append(np.asarray([m0,m12,delta1,delta2,U2sfnuuni,U2sfnumin]))
    xx=np.asarray(xx)

