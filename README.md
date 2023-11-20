This script calculates the problem of a shock followed by a CJ deflagration, 
propagating in a tube with rear venting.

It implements thermal and chemical equilibrium calculated with temperature
dependent heat capacities.

See associated article (to be published in Combustion and Flame, 2023)
for the details on the notation, terminology, method and algorithms.

It is sectionned into 4 different part :

	USER PARAMETERS		Where the user can set the gas initial conditions
						(pressure, temperature, mixture), and the
						rear wall opening area.

	LIST OF FUNCTIONS	Where all the functions used in this script are
						defined.

	PRE-PROCESSING		Where variables derived from the initial state
						are defined.


	MAIN SCRIPT			Where the calculations to determine the flow
						properties are run.


	POST-PROCESSING		Where variables derived from the flow properties
						are calculated and printed.
	

NOTES:

	- This script runs on Python 2. Running it on Python 3 will not work as 
	  the some of the syntax is different.

	- Requires the Python thermo-chemical calculation library CANTERA to work.
	  CANTERA is called when importing library CJ3.
	
	- Requires an older version of the Shock and Detonation Toolbox (SDToolbox)
	  for Python 2 developped by Shephered et al. (2015) to work.
	  It is provided with this package and called in the library CJ3.

	- Library CJ3 contains the functions meant to compute the
	  CJ deflagration state.
	  The method is derived from the CJ detonation calculations
	  implemented in the SDToolbox.


Willstrong Rakotoarison
wrako07a@uottawa.ca
willstrong.rkt@gmail.com

Created on 		Jan. 2019
Last updated on	Nov. 2023
