import pyspheno
import commands
if __name__ == '__main__':
    #Loop over initial x
    #obtain vd and mu for msugra fixed parameters
    pyspheno.LesHouches['SPHENOINPUT'].entries[91]=1
    pyspheno.writeLHAinFile(pyspheno.LesHouches,'LesHouches.in')
    lsout=commands.getoutput(pyspheno.sphenocmd)
    spc,decays=pyspheno.pyslha.readSLHAFile('SPheno.spc')
    chietau=pyspheno.FilterDecay(pyspheno.FilterDecay(decays,pdaug=11),pdaug=15)
    chietaubr=sum([x.br for x in chietau[1000022].decays])
    chimutau=pyspheno.FilterDecay(pyspheno.FilterDecay(decays,pdaug=13),pdaug=15)
    chimutaubr=sum([x.br for x in chimutau[1000022].decays])
    print 'solar mixing %g' %(chietaubr/chimutaubr)

    
#Fron RpParaneters.dat in SPhenoRP_intel/test
#3E-02 0.6 ! epsilon_1
#1.E-05 1. ! epsilon_2
#1.E-05 1. ! epsilon_3
#1.E-05 6.E-02 ! Lambda_1 = v_d epsilon_1 + mu v_L1
#3.E-02 7E-01 ! Lambda_2 = v_d epsilon_2 + mu v_L2
#3.E-02 7E-01 ! Lambda_3 = v_d epsilon_3 + mu v_L3


#
