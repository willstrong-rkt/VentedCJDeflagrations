"""
	This script contains the functions required to calculate the
	CJ deflagration speed from a given initial state, and the 
	final state associated.
	
	It is based on the SDToolbox developped by Shepherd et al.,
	described in the following reference:
	
	Browne, S., Ziegler, J., and Shepherd, J. Numerical Solution Methods for Shock and
	Detonation Jump Conditions, Mar. 2015.
	
	
	Yohan Vilende, Willstrong Rakotoarison
	Contact : wrako027@uottawa.ca
	          willstrong.rkt@gmail.com
	Last modified : November 15, 2023
	
"""


from cantera import *
import cantera as ct
from SDToolbox import *
from scipy.optimize import fmin

# CJ-DEFLAGRATION FUNCTIONS ====================================================

def CJ_calc3(gas, gas1, ERRFT, ERRFV, x, w1_guess=500, T_guess=1000):
    """

    CJ_calc3
    For a reactive discontinuity, this function calculate for a given specific volume 
    (represented by the ratio x = v1/v2), the velocity of the discontinuity and the
    state downstream it.
    This is equivalent to find the point on the Hugoniot curve and the velocity assocuated
    at a given v2 from a known initial state.

    FUNCTION
    SYNTAX
    [gas,w1,T] = CJ_calc3(gas,gas1,ERRFT,ERRFV,x)

    INPUT
    gas = working gas object
    gas1 = gas object at initial state
    ERRFT,ERRFV = error tolerances for iteration
    x = volume ratio(v1/v2)

    OUTPUT
    gas = gas object at equilibrium state
    w1 = initial velocity to yield prescribed density ratio
    T =  temperature of the downstream state
    """
    r1 = gas1.density
    V1 = 1/r1; P1 = gas1.P ; T1 = gas1.T
    i = 0; DT = 1000; DV = 1000; DP = 1000;DW =1000
    
    #PRELIMINARY GUESS
    V = V1/x;
    r = 1/V;
    w1 = w1_guess;
    T = T_guess;
    [P, H] = eq_state(gas,r,T)
    
    #START LOOP
    while (abs(DT) > ERRFT*T and abs(DW) > ERRFV*w1):
        i = i + 1
        if i == 500:
            print 'CJ_calc function could not converge after 500 iterations on the Newton-Raphson algorithm.'
            print "Results given with : abs(DT) - ERRFT*T = %e \t abs(DW) - ERRFW*w1 = %e" \
                   %(abs(DT) - ERRFT*T, abs(DW) - ERRFV*w1)
            break
            #return
            
        #CALCULATE FH & FP FOR GUESS 1
        [FH,FP] = FHFP_CJ(gas,gas1,w1)

        #TEMPERATURE PERTURBATION
        DT = T*0.02; Tper = T + DT;
        Vper = V; Rper = 1/Vper;
        Wper = w1;
        [Pper, Hper] = eq_state(gas,Rper,Tper)
        #CALCULATE FHX & FPX FOR "IO" STATE
        [FHX,FPX] = FHFP_CJ(gas,gas1,Wper)
        #ELEMENTS OF JACOBIAN
        DFHDT = (FHX-FH)/DT; DFPDT = (FPX-FP)/DT;

        #VELOCITY PERTURBATION
        DW = 0.02*w1; Wper = w1 + DW;
        Tper = T; Rper = 1/V;
        [Pper, Hper] = eq_state(gas,Rper,Tper)
        #CALCULATE FHX & FPX FOR "IO" STATE
        [FHX,FPX] = FHFP_CJ(gas,gas1,Wper)
        #ELEMENTS OF JACOBIAN
        DFHDW = (FHX-FH)/DW; DFPDW = (FPX-FP)/DW;

        #INVERT MATRIX
        J = DFHDT*DFPDW - DFPDT*DFHDW
        b = [DFPDW, -DFHDW, -DFPDT, DFHDT]
        a = [-FH, -FP]
        DT = (b[0]*a[0]+b[1]*a[1])/J; DW = (b[2]*a[0]+b[3]*a[1])/J;
    
        #CHECK & LIMIT CHANGE VALUES
        #VOLUME
        DTM = 0.2*T
        if abs(DT) > DTM:
            DT = DTM*DT/abs(DT)
        #MAKE THE CHANGES
        T = T + DT; w1 = w1 + DW;
        [P, H] = eq_state(gas,r,T)

    return [gas, w1, T]
    
def CJspeed3(P1, T1, mix, mech,xmin,xmax,plt_num, ERRFT = 1e-6, ERRFV = 1e-6, w1_guess=500, T_guess=1000):
    """

    CJspeed3
    Calculates CJ deflagration velocity.
    WR: This method is based on the minimum wave speed algorithm

    FUNCTION
    SYNTAX
    [cj_speed,R2,dnew] = CJspeed3(P1,T1,mix,mech,xmin,xmax,plt_num)

    INPUT
    P1 = initial pressure (Pa)
    T1 = initial temperature (K)
    mix = string of reactant species mole fractions
    mech = cti file containing mechanism data (i.e. 'gri30.cti')
    plt_num = unused
    xmin = minimum of volume ratio
    xmax = maximum of volume ratio
    OUTPUT
    cj_speed = CJ deflagration speed (m/s)
    R2 = R-squared value of LSQ curve fit
    dnew = volume ratio at CJ state
    """
    #DECLARATIONS
    numsteps = 20;
    maxv = xmax;
    minv =xmin ;

    w2 = zeros(numsteps,float)
    rr = zeros(numsteps,float)

    gas1 = Solution(mech)
    gas  = Solution(mech)
    
    #INTIAL CONDITIONS
    gas.TPX  = T1, P1, mix
    gas1.TPX = T1, P1, mix     
    
    #INITIALIZE ERROR VALUES & CHANGE VALUES
    #ERRFT = 1*10**-6;  ERRFV = 1*10**-6;

    i = 1;
    T1 = gas1.T;
    P1 = gas1.P;
    
    R2 = 0.0;
    cj_speed = 0.0
    a = 0.0;
    b = 0.0;
    c = 0.0;
    dnew = 0.0

    # Set parameters for the loop
    counter = 1;		# Initialize counter
    mincounter = 4;		# Minimum loop iteration
    maxcounter = 20;	# Maximum loop iteration
    
    # Start loop to find the maximum flame speed
    while (counter <= mincounter) or (R2 < 0.99999):
        step = (maxv-minv)/float(numsteps)
        i = 0
        x = minv
        
		# Calculate several valules of w1 and v1/v2 on the interval [minv;maxv]
        for i in range (0,numsteps) :
            gas.TPX = T1, P1, mix
        
            # Call subfunction to find speed for a specific ratio of volume
            [gas,temp,T] = CJ_calc3(gas, gas1, ERRFT, ERRFV, x, w1_guess, T_guess=T_guess)
            w2[i] = temp 
            rr[i] = gas.density/gas1.density
        
            i = i + 1;
            x = x + step;
		
        # Fit the points of graph w1 = f(v1/v2) to a parabola.
        [a,b,c,R2,SSE,SST] = LSQ_CJspeed(rr,w2)
		
        # Get the locus of maximum(v1/v2) = dnew to be the CJ specific volume.
        dnew = -b/(2.0*a)  
		
        # Define the a new interval where to find maximum flame speed as dnew (+-) 1%.
        minv = dnew - dnew*0.001
        maxv = dnew + dnew*0.001
        
		# Get the current estimation of CJ-deflagration speed being the maximum of value
        # of the parabolic fit.
        cj_speed = a*dnew**2 + b*dnew + c
        
        # Update the counter and re-run the loop until R2 for parabolla fit small enough.		
        counter = counter + 1;
        
        # Stop the loop on maximum speed if too much iterations, and warn user.
        if counter > maxcounter:
            print "\nReached %i iteration on finding maximum flame speed in CJSpeed3 function." %(counter)
            print "Residual on least-square algorithm R2 = ", R2, "\n"
            break
        
    return [cj_speed,R2,dnew]



def CJ_deflagration(P1, T1, mix, mech, ERRFT = 1e-6, ERRFV = 1e-6):
    """

    CJ_deflagration
    Calculates CJ speed for a deflagration , gas at CJ state and Mach Number.

    FUNCTION
    SYNTAX
    [ucj,gas,M] = CJ_deflagration(P1,T1,mix,mech)

    INPUT
    P1 = Initial pressure (Pa)
    T1 = Initial temperature (K)
    mix = string of reactant species mole fractions
    mech = sandiego2014.cti or gri30.cti for example
    
    OPTIONAL ARGUMENTS
    ERRFT = 1e-6 : Tolerance on temperature for the 2-variables (T,v) Newton-Raphson algorithm [K]
    ERRFT = 1e-6 : Tolerance on volume for the 2-variables (T,v) Newton-Raphson algorithm [m3/kg]

    OUTPUT
    ucj = CJ deflagration
    gas = gas object at CJ state
    M =  Mach Number
    """
    
    # Initialize initial gas at  state T1, P1, mix (input)
    gas1 = ct.Solution(mech)
    gas1.TPX = T1,P1,mix

    # Estimate CJ-deflagration volume for perfect gas, and volume ratio associated.
    # This is to be used to find the minimum volume ratio of the interval where to find
    # the CJ-deflagration volume.
    [gamma1, rgas, Qgas] = perfect_gas_properties(ct.one_atm, 298.15, mix, mech)
    [Dcj_id, Pcj_id, Tcj_id, vcj_id, mflowcj_id] = CJ_deflagration_perfect(P1, T1, gamma1, rgas, Qgas)
    
    rho1_id = P1 / (rgas*T1)
    v1_id = 1/rho1_id
    xcj_id = v1_id / vcj_id
    
    # Define the minimum volume ratio of the interval where to find the CJ-deflagration
    # volume, as (xcj_id - 50% of xcj_id), equivalent to take vmax = 2*vcj_id
    xmin = xcj_id * (1-0.5)
    
    # Calculate the constant pressure combustion, that will give the volume of the 
    # post-constant pressure combustion state.
    gas_CP = ct.Solution(mech)
    gas_CP.TPX = T1,P1,mix    
    gas_CP.equilibrate('HP')
    
	# Define the maximum volume ratio of the interval where to find the CJ-deflagration
    # volume, as the one given for CP-combustion volume (minus 1% for safety).
    xmax = (gas1.v / gas_CP.v) * (1-0.01)
    
    # Call function CJspeed3 to find CJ deflagration velocity and x_cj the
    # volume ratio at the CJ deflagration state
    # Uses w1_guess = Dcj_id + 50% to favorise convergence from the top
    [D_cj, R2, x_cj] = CJspeed3(P1, T1, mix, mech, xmin, xmax, 0, ERRFT, ERRFV, w1_guess=Dcj_id*1.5, T_guess=Tcj_id*0.8)
    
    # Calculate gas_cj a gas object at CJ deflagration state
    gas_cj  = ct.Solution(mech)
    gas_cj.TPX  = T1,P1,mix
    [gas_cj, temp, T_cj] = CJ_calc3(gas_cj, gas1, ERRFT, ERRFV, x_cj, w1_guess=Dcj_id*1.5, T_guess=Tcj_id*0.8)
    
    # Calculate the equilibrium sound speed in the CJ-deflagration state to get the 
    # Mach number of the flow in the burnt gases
    [c_cj_eq, c_cj_fr] = equilSoundSpeeds(gas_cj)
    
    # Calculate flow speed and Mach number in the CJ-deflagration state
    uflow_cj = D_cj / x_cj
    Mflow_cj = uflow_cj / c_cj_eq

    return [D_cj, gas_cj, Mflow_cj]







# PERFECT GAS FUNCTIONS =============================================================

def Mcj_def(Pi, vi, Q, gamma):
    """
    
	Mcj_def
    Calculates the CJ deflagration Mach number in the perfect gas approximation
    
    SYNTAX
    Mcj = Mcj_def(Pi, vi, Q, gamma) 
    
    INPUT
    Pi : Initial pressure [Pa]
    vi : Initial volume [m3]
    Q : Energy of combustion of the mixture [J/kg of mixture], estimated for a perfect gas 
    gamma : specific heat ratio of the gas [-], estimated for a perfect gas.
    
    OUTPUT
    Mcj : Deflagration Mach number estimated for a perfect gas
    
    NOTES
    - gamma can be estimated using function "perfect_gas_properties" that returns
      the perfect gas properties of a mixture at a given state.
    
    """
    n = (gamma**2.0-1)/gamma * (Q / (Pi*vi))
    return ( 1 + n - ( (n + 1)**2.0 - 1 )**0.5 )**0.5
    
    

def perfect_gas_properties(P0, T0, mix, mech):
    """ 
    
	perfect_gas_properties
    Uses Cantera to estimate perfect gas properties (gamma, rgas, Qcomb)
    at initial state defined by P0, T0, mix.
    
    SYNTAX
    [gamma, rgas, Qgas] = perfect_gas_properties(P0, T0, mix, mech)
    
    INPUT
    P0 = Initial pressure (Pa)
    T0 = Initial temperature (K)
    mix = string of reactant species mole fractions
    mech= sandiego2014.cti or gri30.cti for example

    OUTPUT
    gamma = Heat capacity ratio
    rgas = specific perfect gas constant
    Qgas = energy of reaction for the gas
    
    """
    
	# Set gas at standard state
    gas = ct.Solution(mech)
    gas.TPX = T0, P0, mix

    # Extract fresh gases properties
    rgas = ct.gas_constant / gas.mean_molecular_weight
    gamma = gas.cp_mass/gas.cv_mass
    uf0 = gas.int_energy_mass

    # Equilibrate at UV constant and bring to T0, P0
    gas.equilibrate("UV")
    gas.TP = T0, P0
    ub0 = gas.int_energy_mass

    # Calculate heat of reaction per kg of fuel
    Qgas = (uf0 - ub0)

    return [gamma, rgas, Qgas]
    
    
def CJ_deflagration_perfect(P1, T1, gamma1, rgas, Qgas):
    """ 

    CJ_deflagration_perfect
    Calculates the CJ_deflagration speed in the approximation of a perfect gas.
    Uses Cantera to estimate perfect gas properties (gamma, Qcomb, rgas).
    
    SYNTAX
    [Dcj, Pcj, Tcj, vcj, Mflow_cj] = CJ_deflagration_perfect(P1, T1, gamma1, rgas, Qgas)
    
    INPUT
    P1 = Initial pressure (Pa)
    T1 = Initial temperature (K)
    gamma1 = Specific heat ratio
    rgas = Perfect gas constant of the gas
    Qgas = Energy of reaction ofhte gas

    OUTPUT
    Dcj = CJ deflagration velocity
    Pcj = post CJ deflagration pressure
    Tcj = post CJ deflagration temperature
    vcj = post CJ deflagration specific volume
    Mflow_cj =  Mach Number of the flow at the post CJ deflagration state
    
    """
    
    # Define initial state
    rho1 = P1 / (rgas*T1)
    v1 = 1/rho1
    c1 = (gamma1 * rgas * T1)**0.5

    # Calculate CJ deflagration speed
    Mcj = Mcj_def(P1, v1, Qgas, gamma1)

    # Calculate vcj and Pcj
    vcj = v1 * (1 + gamma1 * Mcj**2) / (Mcj**2 * (1+gamma1))
    Pcj = P1 * (1 + gamma1 * Mcj**2) / (1 + gamma1)
	
    # Determine the rest of variables
    Tcj = Pcj * vcj / rgas
    ccj = (gamma1 * rgas * Tcj)**0.5
	
    Dcj = Mcj * c1
    uflow_cj = Dcj * vcj/v1
    Mflow_cj = uflow_cj / ccj
	
    return [Dcj, Pcj, Tcj, vcj, Mflow_cj]
