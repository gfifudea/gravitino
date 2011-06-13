#!/usr/bin/env python
import numpy as np


#call spheno with fit on to obtain mu and vd

import commands
lsout=commands.getoutput('ls')

vd=200.

mu=400.
sgn=(-1)**np.random.random_integers(1,2,3)
eps=np.random.uniform(-1,1,3)
Lambda=np.random.uniform(-1,1,3)
vi=(Lambda-eps*vd)/mu

#3E-02 0.6 ! epsilon_1
#1.E-05 1. ! epsilon_2
#1.E-05 1. ! epsilon_3
#1.E-05 6.E-02 ! Lambda_1 = v_d epsilon_1 + mu v_L1
#3.E-02 7E-01 ! Lambda_2 = v_d epsilon_2 + mu v_L2
#3.E-02 7E-01 ! Lambda_3 = v_d epsilon_3 + mu v_L3


#
