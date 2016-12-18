
############################################################
# Parameters file for the main script
# @bayronportilla-2016


############################################################
# Bulk properties of the system

name     = 'Gliese 581 d'  # Name of the planet
M_s      = 0.31            # Stellar mass, in solar masses 
M_p      = 7.10            # Planetary mass, in earth masses 
R        = 1.70            # Planetary radius, in earth radius
BmAC     = 5.0e-5          # Triaxiality, dimensionless
rigidity = 8.0e10          # Unrelaxed rigidity, in Pascals
tau      = 50              # Maxwell (Andrade) time, in years
alpha    = 2               # Andrade's exponent, dimensionless  
e        = 0.27            # Orbital eccentricity, dimensionless
a        = 0.218           # Semimajor axis, in astronomical units
E0       = 0.0             # Initial eccentric anomaly, in degress



############################################################
# Integration parameters

t_ini    = 0.0            # Starting time for simulation, in years
t_end    = 9.00           # Ending time of the simulation, in years
N        = 5000           # (Default) Number of lines to write in the output file 



############################################################
# Initial conditions

theta_ini = 0.0               # Initial sidereal angle, in degrees
p         = 2.51              # Initial resonance order, p = Omega_ini/n                     
a_ini     = a                 # Initial semimajor axis, in astronomical units
e_ini     = e                 # Initial orbital eccentricity, dimensionless
