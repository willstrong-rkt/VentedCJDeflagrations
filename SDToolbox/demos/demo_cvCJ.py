'''
 Shock and Detonation Toolbox
 http://www.galcit.caltech.edu/EDL/public/sdt/SD_Toolbox/

 Use the cv_CJ function to make output files
 for a constant volume explosion simulation where the
 shock front travels at the CJ speed
 
 FIRST - SPECIFY INITIAL CONDITIONS
# P1 = Initial Pressure
# T1 = Initial Temperature
# U1 = Shock Speed (only for speficied shock speed calculation, not needed
#         for the CJ case)
# q = Initial Composition MUST use capital leters - STRING
# mech = Mechanism File name is CTI format - STRING
#         (Generally all mechanism files are stored in
#         'CANTERA/data' on a Linux machine where CANTERA is the
#         cantera path
#
# NOW - RUN PROGRAM OF CHOICE
# ARGUMENTS
#  ARG 1: Figure Number for CV profile
#  ARG 2: Initial Pressure -- P1
#  ARG 3: Initial Temperature -- T1
#  ARG 4: Initial Composition -- q
#  ARG 5: Mechanism File -- mech
#  ARG 6: Output file name - Enter a string to generate a .plt file
#         starting with that string OR enter '0' (the number zero) to
#         bypass generating an output file
'''

# Updated 2-20-2014 for compatibility with
# Python 2.7 and Cantera 2.1.
from SDToolbox import *

P1 = 100000; T1 = 300; U1 = 2000; P1atm = P1/one_atm;
q = 'H2:2 O2:1 N2:3.76'
mech = 'h2air_highT.cti'

[cj_speed,gas] = cv_CJ(0, P1, T1, q, mech, 'h2air')

print ' '
print 'For ' + q + ' with P1 = %.2f atm & T1 = %.2f K using ' % (P1atm,T1) + mech 
print 'CJ Speed is %.2f m/s' % cj_speed
print 'Final State (after Constant Volume Explosion)'
gas()
print ' '

print 'DONE'
