#!/usr/bin/env python
'''
* ranlog: for generate random numbers in ranges of several orders
  of magnitude
* Functions to analyse neutrino solutions
'''
from pyspheno import * #includes pyslha
import numpy as np
import commands
import sys
import matplotlib.pyplot as plt

def oscilation(spcfile):
    """oscilation parameters
    spcfile is used instead of the dictionary in order to be used
    independtly. This sould delay a little bit the main program"""
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
    spcfile is used instead of the dictionary in order to be used
    independtly. This sould delay a little bit the main program
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

def chisq(spcfile):
    (Delta2m32,Delta2m21,s223,s212,U13)=oscilation(spcfile)
    #constants going into the chisq:
    s212r=0.304;s223r=0.5
    U213=0.001
    m22=2.40;m21=7.65
    if U13**2<0.056:
        return ((s212-s212r)/0.07)**2+((s223-s223r)/0.165)**2\
           +((1E3*Delta2m32-m22)/0.4)**2+((1E5*Delta2m21-m21)/0.6)**2
    else:
        return 10E6
    
def contraints(spcfile):
    '''spcfile is used instead of the dictionary in order to be used
    independtly. This sould delay a little bit the main program'''
    spc,decays=pyslha.readSLHAFile(spcfile)
    check=False
    chi_20=spc['MASS'].entries[1000023]
    mh=spc['MASS'].entries[25]
    bsgamma=spc['SPHENOLOWENERGY'].entries[1]
    deltagmu_2=2.*spc['SPHENOLOWENERGY'].entries[21]
    brmupmum=spc['SPHENOLOWENERGY'].entries[5]
    if chi_20>103.5 and mh>93.5 and (bsgamma>2.77E-4 and bsgamma<4.37E-4)\
           and (deltagmu_2>-11.4E-10 and deltagmu_2<9.4E-9)\
           and brmupmum<4.2E-8:
        check=True
        
    return check

def smparameters(spcfile):
    '''Calculation of s2w from standard model relations:
    See 0705.4263 (4.15)'''
    spc,decays=pyslha.readSLHAFile(spcfile)
    MW=spc['MASS'].entries[24]
    alphaeinv=spc['SMINPUTS'].entries[1]
    GF=spc['SMINPUTS'].entries[2]
    v=1./np.sqrt(np.sqrt(2)*GF)
    s2w=np.pi/(alphaeinv*MW**2*np.sqrt(2)*GF)
    g2=2*MW/v
    return s2w,g2
    
def analyticalU2sgnu(spcfile):
    '''g2 calculated from 0705.4263
    spcfile is used instead of the dictionary in order to be used
    independtly. This sould delay a little bit the main program'''
    spc,decays=pyslha.readSLHAFile(spcfile)
    M=np.asarray([spc['EXTPAR'].entries[i] for i in range(4)])
    nvd=spc['SPHENORP'].entries[15]
    nvu=spc['SPHENORP'].entries[16]
    lamb=np.asarray(spc['SPHENORP'].entries.values()[0:3])
    epsi=np.asarray(spc['RVKAPPA'].entries.values())
    vi=np.asarray(spc['RVSNVEV'].entries.values())
    nmu=((lamb-nvd*epsi)/vi)[0]
    s2w,g2=smparameters(spcfile)
    sw=np.sqrt(s2w)
    cw=np.sqrt(1-s2w)
    g1=g2*sw/cw
    chimatrix=[[M[1],0,-g1*nvd/2.,g1*nvu/2.],[0,M[2],g2*nvd/2.,-g2*nvu/2.],\
               [-g1*nvd/2.,g2*nvd/2,0,-nmu],[g1*nvu/2,-g2*nvu/2,-nmu,0]]
    return nmu**2*s2w*g2**2*(M[2]-M[1])**2/(4*np.linalg.det(chimatrix)**2)*np.linalg.norm(lamb)**2
    
def plotU2sgnu_alpha(x,dots='ro',anal=False):
    x=np.asarray(x)
    if anal:
        #too keep compatibility with old data: :,-1
        return plt.loglog(np.abs(x[:,2]/x[:,1]),x[:,-1],dots)
    else:
        return plt.loglog(np.abs(x[:,2]/x[:,1]),x[:,4],dots)

def finalplot(frame=False):
    plt.xlim(1E-1,5E2)
    plt.ylim(1E-17,2E-11)
    w=np.load('universal.npy')
    x=np.load('nonuniversal12.npy')
    y=x[np.abs(x[:,7])>50]
    zz=np.load('nonuniversaldelta2_0.npy')
    z=zz[np.logical_and(zz[:,5]<1,zz[:,5]>-1)]
    if frame:
        #just the frame and universal point.
        #The basis to draw the regions in inkscape.
        plotU2sgnu_alpha(w,dots='r*')
    else:
        plotU2sgnu_alpha(x,dots='ro')
        plotU2sgnu_alpha(y,dots='yo')
        plotU2sgnu_alpha(z,dots='ko')
        
    plt.hlines(w[:,4],plt.xlim()[0],w[:,2]/w[:,1],linestyles='--')
    plt.vlines(w[:,2]/w[:,1],plt.ylim()[0],w[:,4],linestyles='--')
    plt.xlabel(r'$M_2/M_1$',size=20)
    plt.ylabel(r'$|U_{\tilde{\gamma}\nu}|^2$',size=20)

    plt.show()


def zoom(frame=False):
    x=np.load('nonuniversal12.npy')
    zz=np.load('nonuniversaldelta2_0.npy')
    z=zz[np.logical_and(zz[:,5]<1,zz[:,5]>-1)]
    if not frame:
        plotU2sgnu_alpha(x,dots='ro')
        plotU2sgnu_alpha(z,dots='ko')
    else:
        plotU2sgnu_alpha([z[0]],dots='ko')

    plt.xscale('linear')
    plt.xlim(0.75,1.2)
    plt.ylim(1E-17,1E-14)
    plt.xlabel('')
    plt.ylabel('')

def U2sgnuf(x):
    Ni1=x[0:3]
    Ni2=x[3:6]
    thw=x[6]
    return np.dot(np.cos(thw)*Ni1+np.sin(thw)*Ni2,\
                              np.cos(thw)*Ni1+np.sin(thw)*Ni2)
    
if __name__ == '__main__':
    '''Variables
    xx[:,0]=Q; xx[:,1]=M1, xx[:,2]=M2, xx[:,3]=M3
    xx[:,4]=U2sgnu,xx[:,5]=delta1,xx[:,6]=delta2,
    xx[:,7]=chi_10,xx[:,8]=chi_20,xx[:,9]=chisq
    x[:,10]=U2sgnu_analytical'''
    spcfile='SPheno.spc'
    #sphenocmd='SPheno_intel'
    sphenocmd='SPheno3.1.4' 
    ndelta1=1 #100 #30
    ndelta2=1 #100
    delta1min=0; delta1max=0 #-1,3 for delta2=0
    delta2min=0; delta2max=0 #-1,1
    xx=[]
    m0=1000.
    m12=500.
    #Check universal point and get SM parameters:
    LesHouches=buildLHAinFile()
    LesHouches['MINPAR'].entries[1]=m0
    LesHouches['MINPAR'].entries[2]=m12
    writeLHAinFile(LesHouches,'LesHouches.in')
    lsout=commands.getoutput(sphenocmd)
    if not lsout:
        s2w,g2=smparameters(spcfile)
        thw=np.arcsin(np.sqrt(s2w))
    else:
        print "Wrong SUGRA point"
        sys.exit()

    #prepare non-sugra point
    LesHouches=buildLHAinFile(universal=False) #see below
    j=0
    for delta1 in np.random.uniform(delta1min,delta1max,ndelta1):
        for delta2 in np.random.uniform(delta2min,delta2max,ndelta2):
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
                M=np.asarray([spc['EXTPAR'].entries[i] for i in range(4)])
                Ni1=np.asarray([spc['RVNMIX'].entries[i][1] for i in range(1,4)])
                Ni2=np.asarray([spc['RVNMIX'].entries[i][2] for i in range(1,4)])
                xin=np.hstack((Ni1,Ni2,thw))
                U2sgnu=U2sgnuf(xin)
                chi_10=spc['MASS'].entries[1000022]
                chi_20=spc['MASS'].entries[1000023]
                appU2sgnu=analyticalU2sgnu(spcfile)
                other=np.asarray([U2sgnu,delta1,delta2,chi_10,chi_20,\
                                  chisq(spcfile),appU2sgnu])
                if contraints(spcfile):
                    xx.append(np.hstack((M,other)))
                else:
                    if chi_20>103.5: print 'const',delta1,delta2

            else:
                print 'bad',delta1,delta2

            if j%1000==0: print j
            j=j+1

    xx=np.asarray(xx)
    np.save('nonuniversal',xx)
    #xx=np.load('nonuniversal.npy')
    #plt.loglog(np.abs(xx[:,2]/xx[:,1]),xx[:,4],'ro')
    #plotU2sgnu_alpha(xx)
    plt.show()
