'''
 Shock and Detonation Toolbox
 http://www.galcit.caltech.edu/EDL/public/sdt/SD_Toolbox/

This demo generates data in output files to plot the Rayliegh and Hugoniot
lines.

These output files are in the tecplot format (the header specifically),
however they can also be easily read into excel for plotting


Creates arrays for Rayleigh Line with slope U1, Reactant,
and Product Hugoniot Curves

# FIRST - SPECIFY INITIAL CONDITIONS & DESIRED PLOTS
# P1 = Initial Pressure
# T1 = Initial Temperature
# q = Initial Composition MUST use capital letters - STRING
# mech = Mechanism File name in CTI format - STRING
#       (Generally all mechanism files are stored in 'Cantera/data'
#        on a Linux Machine)
# fname = Mixture name used for output file names (use zero '0' if
#         no output file desired
# od = Overdrive for Rayleigh Line -- multiplied by CJ speed to find the
#      slope of the Rayleigh Line
# gbound = Overdrive for Gamma Bound -- mulitiplied by CJ speed to find
#          the maximum gamma
#
'''

# Updated 2-20-2014 for compatibility with
# Python 2.7 and Cantera 2.1.

from SDToolbox import *
from scipy.optimize import fsolve

P1 = 100000; T1 = 300;
q = 'H2:2 O2:1 N2:3.76'
mech  = 'h2air_highT.cti'
fname = 'h2air'
od = 1; gbound = 1.4

# EDIT VALUES ABOVE THIS LINE
##############################
gas1 = Solution(mech);
gas1.TPX = T1, P1, q


h1 = gas1.enthalpy_mass
r1 = gas1.density
v1 = 1.0/gas1.density

#Get CJ Point
[cj_speed,R2] = CJspeed(P1,T1,q,mech,0)
gas = PostShock_eq(cj_speed,P1,T1,q,mech)
vcj = 1.0/gas.density
Pcj = gas.P/101325.0

print 'CJ Point Found'

U1 = od*cj_speed

#Find Postshock specific volume for U1
gas = PostShock_fr(U1, P1, T1, q, mech)
vsj = 1.0/gas.density
Psj = gas.P/101325.0

#Find Gamma 
gas = PostShock_fr(gbound*cj_speed, P1, T1, q, mech)
g = gas.cp_mass/gas.cv_mass

# RAYLEIGH LINE, SJUMP MIN & MAX, GAMMA CONSTRAINT
minp = 0.1*vcj; maxp = 2.0*vcj; step = 0.01*vcj;
i = 0; v2 = minp; 

##### NUMPY?
n = int((maxp-minp)/step);
vR = zeros(n,float);
PR = zeros(n,float); 

min_line = zeros(n,float);
max_line = zeros(n,float);
gamma = zeros(n,float); 

while(i < n):
    vR[i] = v2;
    PRpa = (P1 - r1**2*U1**2*(v2-v1))
    PR[i] = PRpa/101325.0

    min_line[i] = 1.0/40.0; max_line[i] = 1.0/1.005;
    gamma[i] = (g-1.0)*r1/(g+1.0)
    
    i = i + 1; v2 = v2 + step
print 'Rayleigh Line Array Created'

# REACTANT HUGONIOT
gas.TPX = T1, P1, q; 
Ta = gas.T; va = 1.0/gas.density

n = int((va-vsj)/0.01);
PH1 = zeros(n+1,float);
vH1 = zeros(n+1,float);

PH1[0] = gas.P/101325.0; vH1[0] = va

i = 0; vb = va; 
while(i < n):
    vb = va - (i+1)*0.01;
    #Always starts at new volume and previous temperature to find next state
    fval = fsolve(hug_fr,Ta,args=(vb,h1,P1,v1,gas))
    gas.TD = fval, 1.0/vb
    PH1[i+1] = gas.P/101325.0; vH1[i+1] = vb
    i = i + 1

print 'Reactant Hugoniot Array Created'

gas = PostShock_eq(U1, P1, T1, q, mech)
veq = 1.0/gas.density
Peq = gas.P/101325.0

# PRODUCT HUGONIOT

# Get the first point on the product Hugoniot - CV combustion
gas.TPX = T1, P1, q
gas.equilibrate('UV')
Ta = gas.T
va = 1.0/gas.density

n = int((va-0.4*vcj)/0.01);
PH2 = zeros(n+1,float);
vH2 = zeros(n+1,float);

PH2[0] = gas.P/101325.0; vH2[0] = va;

i = 0; vb = va;
while(i < n):
    vb = va - (i+1)*0.01
    #Always starts at new volume and previous temperature to find next state
    fval = fsolve(hug_eq,Ta,args=(vb,h1,P1,v1,gas))
    gas.TD = fval, 1.0/vb
    PH2[i+1] = gas.P/101325.0; vH2[i+1] = vb
    i = i + 1

print 'Product Hugoniot Array Created'

if(fname == 0):
    print 'No output files created (RLineHug_script)'
else:
    #Create Output Data File for Tecplot
    import datetime
    a = size(PR)-1; k = 0;
    d = datetime.date.today(); P = P1/101325.0

    fn = fname + '_%d_RgLines.plt' % U1
    outfile = file(fn, 'w')
    outfile.write('# RAYLEIGH LINE, SJUMP MAX & MIN, AND GAMMA CONSTRAINT FOR GIVEN SHOCK SPEED\n')
    outfile.write('# CALCULATION RUN ON %s\n\n' % d)
    outfile.write('# INITIAL CONDITIONS\n')
    outfile.write('# TEMPERATURE (K) %.1f\n' % T1)
    outfile.write('# PRESSURE (ATM) %.1f\n' % P)
    outfile.write('# DENSITY (KG/M^3) %.4f\n' % r1)
    outfile.write('# SPECIES MOLE FRACTIONS: ' + q + '\n')
    outfile.write('# MECHANISM: ' + mech + '\n')
    outfile.write('# SHOCK SPEED (M/S) %.2f\n\n' % U1)
    outfile.write('# THE OUTPUT DATA COLUMNS ARE:\n')
    outfile.write('Variables = "Specific Volume", "Pressure", "Max Line", "Min Line", "Gamma Line"\n')
    while k <= a:
        outfile.write('%.4E \t %.1f \t %.4f \t %.4f \t %.4f\n' % (vR[k], PR[k], min_line[k], max_line[k], gamma[k]))
        k = k + 1
    
    a = size(PH1)-1; k = 0;
    fn = fname + '_ReacHug.plt' 
    outfile = file(fn, 'w')
    outfile.write('# REACTANT HUGONIOT FOR GIVEN MIXTURE\n')
    outfile.write('# CALCULATION RUN ON %s\n\n' % d)
    outfile.write('# INITIAL CONDITIONS\n')
    outfile.write('# TEMPERATURE (K) %.1f\n' % T1)
    outfile.write('# PRESSURE (ATM) %.1f\n' % P)
    outfile.write('# DENSITY (KG/M^3) %.4f\n' % r1)
    outfile.write('# SPECIES MOLE FRACTIONS: ' + q + '\n')
    outfile.write('# MECHANISM: ' + mech + '\n\n')
    outfile.write('# THE OUTPUT DATA COLUMNS ARE:\n')
    outfile.write('Variables = "Specific Volume", "Pressure"\n')
    while k <= a:
        outfile.write('%.4E \t %.1f\n' % (vH1[k], PH1[k]))
        k = k + 1
    
    a = size(PH2)-1; k = 0;
    fn = fname + '_ProdHug.plt' 
    outfile = file(fn, 'w')
    outfile.write('# PRODUCT HUGONIOT FOR GIVEN MIXTURE\n')
    outfile.write('# CALCULATION RUN ON %s\n\n' % d)
    outfile.write('# INITIAL CONDITIONS\n')
    outfile.write('# TEMPERATURE (K) %.1f\n' % T1)
    outfile.write('# PRESSURE (ATM) %.1f\n' % P)
    outfile.write('# DENSITY (KG/M^3) %.4f\n' % r1)
    outfile.write('# SPECIES MOLE FRACTIONS: ' + q + '\n')
    outfile.write('# MECHANISM: ' + mech + '\n\n')
    outfile.write('# THE OUTPUT DATA COLUMNS ARE:\n')
    outfile.write('Variables = "Specific Volume", "Pressure"\n')
    while k <= a:
        outfile.write('%.4E \t %.1f\n' % (vH2[k], PH2[k]))
        k = k + 1    
    
    fn = fname + '_%d_RLineHug_points.plt' % U1 
    outfile = file(fn, 'w')
    outfile.write('# SPECIFIC POINTS FOR GIVEN SHOCK SPEED\n')
    outfile.write('# CALCULATION RUN ON %s\n\n' % d)
    outfile.write('# INITIAL CONDITIONS\n')
    outfile.write('# TEMPERATURE (K) %.1f\n' % T1)
    outfile.write('# PRESSURE (ATM) %.1f\n' % P)
    outfile.write('# DENSITY (KG/M^3) %.4f\n' % r1)
    outfile.write('# SPECIES MOLE FRACTIONS: ' + q + '\n')
    outfile.write('# MECHANISM: ' + mech + '\n')
    outfile.write('# SHOCK SPEED (M/S) %.2f\n\n' % U1)
    outfile.write('# THE OUTPUT DATA COLUMNS ARE:\n')
    outfile.write('# ORDER = Initial, PostShock, CJ, Equilibrium\n')
    outfile.write('Variables = "Specific Volume", "Pressure"\n')
    outfile.write('%.4E \t %.1f\n' % (v1, P))
    outfile.write('%.4E \t %.1f\n' % (vsj, Psj))
    outfile.write('%.4E \t %.1f\n' % (vcj, Pcj))
    outfile.write('%.4E \t %.1f\n' % (veq, Peq))
