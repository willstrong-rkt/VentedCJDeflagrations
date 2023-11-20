'''
 Shock and Detonation Toolboox^
 http://www.galcit.caltech.edu/EDL/public/sdt/SD_Toolbox/

Calculate the CJ speed using the Minimum Wave Speed Method

This CJspeed function is more robust (works with large mechanisms) but slower
than the CJspeed2 function (based on the equilibrium sound speed algorithm)

Inputing a 1 (rather than zero) to the CJspeed function generates plots
that show the goodness fit of the minimization

'''

# Updated 2-20-2014 for compatibility with
# Python 2.7 and Cantera 2.1.

from SDToolbox import *

P1 = 100000; P1atm = P1/one_atm;
T1 = 300;
U = 2000;
q = 'H2:2 O2:1 N2:3.76'
mech = 'h2air_highT.cti'
 
[cj_speed,R2] = CJspeed(P1, T1, q, mech, 0);   
print ' '
print 'For ' + q + ' with P1 = %.2f atm & T1 = %.2f K using ' % (P1atm,T1) + mech 
print 'CJ Speed is %.2f m/s' % cj_speed
print ' '
