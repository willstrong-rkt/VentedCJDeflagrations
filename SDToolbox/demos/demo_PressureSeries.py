'''
 Shock and Detonation Toolbox
 http://www.galcit.caltech.edu/EDL/public/sdt/SD_Toolbox/

 This is a demostration of how to vary the initial Pressure in a loop for
 constant volume explosions 

 Here the explosion functions are called for varying conditions,
 new properties are calculated, plots are displayed, and output files are
 written

 Using this demo as a guide, users can design their own loops and calulations
'''
from SDToolbox import *
print 'Initial Pressure Series'

mech = 'h2air_highT.cti';   
gas  = Solution(mech); 
gas1 = Solution(mech);
nsp  = gas.n_species;

npoints=15;
P1          = zeros(npoints,float)
cj_speed    = zeros(npoints,float)
exo_time_CV = zeros(npoints,float)
ind_time_CV = zeros(npoints,float)
theta_effective_CV = zeros(npoints,float)
T2   = zeros(npoints,float)
P2   = zeros(npoints,float)
rho2 = zeros(npoints,float)
Ts   = zeros(npoints,float)
Ps   = zeros(npoints,float)

# find Hydrogen nitrogen, and oxygen indices
ih2  = gas.species_index('H2');
io2  = gas.species_index('O2');
in2  = gas.species_index('N2');
x = 'H2:2.0,O2:1.0,N2:3.76';   
T1 = 300
Po = 100000;
fig_num = 0;
fname = 0;

print 'Initial Conditions'
print x + ', Temperature = %.2f K' % (T1)
print "For %s initial pressures" % (npoints)
i = 0;
while i  < npoints:
    P1[i] = Po*(0.1 +1.1/npoints*(i));
    P = P1[i];
    print '%i : P1 = %.2f atm' % (i,P/one_atm)

    gas.TPX = T1,P1[i],x;
    
    ###Constant Volume Explosion Data###
    # FIND POST SHOCK STATE FOR GIVEN SPEED
    [gas,cj_speed[i]] = CJspeed2(P1[i], T1, x, mech);
    gas = PostShock_fr(cj_speed[i], P1[i], T1, x, mech);
    Ts[i] = gas.T; #frozen shock temperature   
    Ps[i] = gas.P; #frozen shock pressure
    # SOLVE CONSTANT VOLUME EXPLOSION ODES
    b = 10000; j = gas.n_species;
    out = cvoutput(b,j)
    out = explosion(gas,fname,out);
    exo_time_CV[i] = out.exo_time;
    ind_time_CV[i] = out.ind_time;
    ##Calculate CJstate Properties###
    gas = PostShock_eq(cj_speed[i],P1[i], T1,x, mech);
    T2[i] = gas.T;
    P2[i] = gas.P;
    rho2[i] = gas.density;
   
    #Approximate the effective activation energy using finite differences
    Ta = Ts[i]*(1.01);
    gas.TPX = Ta,Ps[i],x; 
    out = explosion(gas,0,out);
    taua = out.ind_time;

    Tb = Ts[i]*(0.99);
    gas.TPX = Tb,Ps[i],x; 
    out = explosion(gas,0,out);
    taub = out.ind_time;
    #Approximate effective activation energy for CV explosion
    theta_effective_CV[i] = 1.0/Ts[i]*((log(taua)-log(taub))/((1.0/Ta)-(1.0/Tb)));   
    i = i + 1;

##################################################################################################
# CREATE OUTPUT TEXT FILE
##################################################################################################
d = datetime.date.today(); 
fn = 'Pressure_Series.plt';
outfile = file(fn, 'w');
outfile.write('# CONSTANT VOLUME EXPLOSION\n');
outfile.write('# CALCULATION RUN ON %s\n\n' % d);
outfile.write('# INITIAL CONDITIONS\n');
outfile.write('# TEMPERATURE (K) %0.2f\n' % T1);
q = 'H2:2 O2:1 N2:3.76';
outfile.write( '# SPECIES MOLE FRACTIONS: %s \n' % q );
outfile.write('# THE OUTPUT DATA COLUMNS ARE:\n');
outfile.write('Variables = "Pressure state 1", "Temperature state 2", "Pressure state 2", "density state 2", "CJ Speed", "Temperature Post Shock", "Pressure Post Shock", "Induction time CV", "Exothermic time CV",  "Effective Activation Energy CV "\n');
for i in range(npoints):
     outfile.write('%.4E \t %.4E \t %.4E \t %.4E \t %.4E \t %.4E \t %.4E \t %.4E \t %.4E \t %.4E \n'% (P1[i], T2[i], P2[i], rho2[i], cj_speed[i], Ts[i], Ps[i], ind_time_CV[i], exo_time_CV[i], theta_effective_CV[i]));


