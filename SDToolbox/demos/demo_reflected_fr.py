"""
 Shock and Detonation Toolbox
 http://www.galcit.caltech.edu/EDL/public/sdt/SD_Toolbox/

Demo Reflected Shock - Calculates Relected Post-shock state for a specified shock speed
 speed for a specified mixture (FROZNE not equilibrium calculations)

INITIAL CONDITIONS
P1 = Initial Pressure
T1 = Initial Temperature
UI = Incident Shock Speed 
q = Initial Composition MUST use capital letters - STRING
mech = Mechanism File name in CTI format - STRING 
      (Generally all mechanism files are stored in 'Program Files/Common Files/Cantera/data'
       on a Windows Machine)
OUTPUT
UR = Rflected shock speed for state 1-->3
gas3 = gas object with properties at state 3 
             (between intial state and wall after the reflection)

Demo Reflected Shock Psuedocode
set P1, T1, q, UI
 gas2 = PostShock_fr
 gas3, UR <-- reflected_fr
 call reflected_fr
     u2 = sqrt((p2-p1)*(v1-v2)) ;%particle velocity
     guess P3 for constant gamma
       --> approximate
           UR = (p3-p2)/u2/rho2-u2
           rho3 = (p3-p2)/((p3-p2)/rho2-u2^2)
           T3 = T1*p3*(1/rho3)/(p1*(1/rho1))
      gas3 <-- PostReflectedShock_fr(u2,gas2, gas3)
                while((abs(DT) > ERRFT*T) | (abs(DV) > ERRFV*V))
                    del H, del P  <---  FHFP_reflected_fr(u2,gas3,gas2);
                    perturb the temperature
                    gas3 <--- set state
                    del Hx, del Px  <---  FHFP_reflected_fr
                              FH = H3 -H2- 1/2*u2^2* (r3/r2+1)/(r3/r2-1);
                              FP = P3 - P2 - r3*u2^2/(r3/r2-1);
                    perturb the specific volume
                    gas3 <--- state
                    del Hx, del Px  <---  FHFP_reflected_fr
                    Update T and V
                    gas3 <--- state                 
      UR = (p3-p2)/u2/rho2-u2
"""

# Updated 2-20-2014 for compatibility with
# Python 2.7 and Cantera 2.1.

from SDToolbox import *
import math

P1 = 100000; T1 = 300; UI = 2000; u1 = 0;
P1atm = P1/one_atm;
q = 'H2:2 O2:1 N2:3.76';    
mech = 'h2air_highT.cti';               
print ' '
print 'For ' + q + ' with P1 = %.2f & T1 = %.2f using ' % (P1atm,T1) + mech
print ' '

gas2 = PostShock_fr(UI, P1, T1, q, mech);
P2 = gas2.P/one_atm;
print 'Frozen Post-Incident-Shock State (%.2f m/s)' % (UI)
print 'T2 = %.2f K & P2 = %.2f atm' % (gas2.T,P2)
print ' '

gas1 = Solution(mech);
gas3 = Solution(mech);
gas1.TPX = T1, P1, q

[p3,UR,gas3] = reflected_fr(gas1,gas2,gas3,UI);
P3 = gas3.P/one_atm;
print 'Frozen Post-Reflected-Shock State'
print 'T3 = %.2f K & P3 = %.2f atm' % (gas3.T,P3)
print "Reflected Wave Speed = %.2f m/s" % (UR)
print ' '
