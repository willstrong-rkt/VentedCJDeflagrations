'''
 Shock and Detonation Toolboox
 http://www.galcit.caltech.edu/EDL/public/sdt/SD_Toolbox/

Calculates  Post-shock state and CJ speed for a specified mixture for a CJ detonation using the algorithm based on minumizing the shock speed

NOTE: The function, 'CJspeed', is more robust (works with large mechanisms) but slower than the 'CJspeed2' function (based on the equilibrium sound speed algorithm)

INITIAL CONDITIONS
P1 = Initial Pressure
T1 = Initial Temperature
q = Initial Composition MUST use capital letters - STRING
mech = Mechanism File name in CTI format - STRING
   (Generally all mechanism files are stored in
   'Program Files/Common Files/Cantera/data'   on a Windows Machine)

'''

# Updated 2-20-2014 for compitability with
# Python 2.7 and Cantera 2.1.

from SDToolbox import *


P1 = 100000; P1atm = P1/one_atm;
T1 = 300;
q = 'H2:2 O2:1'
mech = 'h2air_highT.cti'

#Ethylene-Oxygen
P1=one_atm; P1atm=P1/one_atm;
T1=300;
q='C2H4:1 O2:3';
mech = 'gri30_highT.cti';

 
[cj_speed,R2] = CJspeed(P1, T1, q, mech, 0);   

gas = PostShock_eq(cj_speed, P1, T1, q, mech)
Ps = gas.P/one_atm

print ' '
print 'CJ State'
gas()
Ps = gas.P/one_atm

print ' '
print 'For ' + q + ' with P1 = %.2f atm & T1 = %.2f K using ' % (P1atm,T1) + mech 
print 'CJ Speed is %.2f m/s' % cj_speed

print 'The CJ State is %.2f atm & %.2f K' % (Ps,gas.T)

[ae,af] = equilSoundSpeeds(gas)

print 'The sound speeds are: af = %.2f m/s & ae = %.2f m/s' % (af,ae)
print ' '
