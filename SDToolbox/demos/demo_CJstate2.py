"""
 Shock and Detonation Toolbox
 http://www.galcit.caltech.edu/EDL/public/sdt/SD_Toolbox/

Calculates  Post-shock state and CJ speed for a specified mixture for a CJ detonation using the algorithm based one the equilibrium sound speed

NOTE: The other function, 'CJspeed', is more robust (works with large mechanisms) but slower than the 'CJspeed2' function (based on the equilibrium sound speed algorithm)
                
INITIAL CONDITIONS
P1 = Initial Pressure
T1 = Initial Temperature
q = Initial Composition MUST use capital letters - STRING
mech = Mechanism File name in CTI format - STRING 
   (Generally all mechanism files are stored in
   'Program Files/Common Files/Cantera/data'   on a Windows Machine)


--- CJSpeed2 Psuedo Code ---
Initialize gas1
guess gas state at CJ point
Vg = V1/5
equilibrate(gas, 'UV')

while((abs(DT) > ERRFT*T) | (abs(DV) > ERRFV*V))
      del Hx, del Px, cj_speed = FHFP_CJ2
       perturb the temperature
       gas <--- eq_state (constant TV)
       del Hx, del Px, cj_speed  <---  FHFP_CJ2
                   cj_speed = aequil <----equilSoundSpeeds
                                equilibrate(gas, 'TP')
                                get s0, p0, r0 
                                r1 = r0*1.0001 (perturb the density)
                                set(gas, 'S', s0, 'V', 1./r1)
                                get pfrozen
                                 equilibrate(gas,'SV');
                                 get p1
                                 aequil = sqrt((p1 - p0)/(r1 - r0));
                   w2 = aequil;
                   w1^3= w2^2*(r2/r1)^2;
                   FH = H2 + 0.5*w2s - (H1 + 0.5*w1s);
                   FP = P2 + r2*w2s - (P1 + r1*w1s);
        perturb the specific volume
        gas <--- eq_state
        del Hx, del Px, cj_speed  <---   FHFP_CJ2
        Update T and V
        gas <--- eq_state
CJ Speed = w1 (w2 = c2 = UCJ - u2)      
"""

# Updated 2-20-2014 for compatibility with
# Python 2.7 and Cantera 2.1.

from SDToolbox import *

P1 = 100000; P1atm = P1/one_atm;
T1 = 300;
q = 'H2:2 O2:1'
mech = 'h2air_highT.cti'


 
[gas, cj_speed] = CJspeed2(P1, T1, q, mech);

print ' '
print 'CJ State'
gas()
Ps = gas.P/one_atm

print ' '
print 'For ' + q + ' with P1 = %.2f atm & T1 = %.2f K using ' % (P1atm,T1) + mech 
print 'CJ Speed is %.2f m/s' % cj_speed


print 'The CJ State is %.2f atm & %.2f K' % (Ps,gas.T)
print ' '

