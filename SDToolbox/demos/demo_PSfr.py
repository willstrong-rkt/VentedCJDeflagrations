'''
 Shock and Detonation Toolbox
 http://www.galcit.caltech.edu/EDL/public/sdt/SD_Toolbox/

 calculate the post shock state based the initial gas
 properties and the shock speed 

 (This is the FROZEN rather than equilibrium state)

 # FIRST - SPECIFY INITIAL CONDITIONS
# P1 = Initial Pressure
# T1 = Initial Temperature
# U = Shock Speed (only for specified shock calculation, not needed for CJ case)
# q = Initial Composition MUST use capital letters - STRING
# mech = Mechanism File name in CTI format - STRING 
#       (Generally all mechanism files are stored in 'Cantera/data'
#        on a Linux Machine)
#
# NOW - RUN PROGRAM OF CHOICE
# ARGUMENTS
#   ARG 1: Initial Pressure -- P1
#   ARG 2: Initial Temperature -- T1
#   ARG 3: Initial Composition -- q
#   ARG 4: Mechanism File -- mech
#   ARG 5: File Name Initial Velocity as a function of density
#       ratio (0 - if no file requested)
'''

# Updated 2-20-2014 for compatibility with
# Python 2.7 and Cantera 2.1.

from SDToolbox import *

P1 = 100000; P1atm = P1/one_atm;
T1 = 300;
U = 2000;
q = 'H2:2 O2:1 N2:3.76'
mech = 'h2air_highT.cti'
 
gas = PostShock_fr(U, P1, T1, q, mech)
Ps = gas.P/one_atm

print ' '
print 'For ' + q + ' with P1 = %.2f & T1 = %.2f using ' % (P1atm,T1) + mech
print 'Ts = %.2f K & Ps = %.2f atm' % (gas.T,Ps)
print ' '
