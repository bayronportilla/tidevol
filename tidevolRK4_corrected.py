#!/usr/bin/env python

from scipy.integrate import odeint,ode
from scipy.optimize import newton
from scipy.special import binom,gamma
from scipy import arctan,special
import numpy as np
from numpy.linalg import norm
from Units import units
from Ecc_Polynomial import G_function
import time



############################################################
# Constants and definitions

Msol   = 1.988500e30 #kg
Mearth = 5.9723e24 #kg
Rearth = 6371000.0 #m
AU     = 149.6e9 #m

############################################################
# Canonical units

Unidades_Canonicas = units(uM=0.001*Mearth,uL=0.218*AU)
uM = Unidades_Canonicas[0]
uL = Unidades_Canonicas[1]
uT = Unidades_Canonicas[2]



############################################################
# Bulk properties of the system. Each property is given in SI

G          = 6.67e-11             # Gravitational constant 
name       = 'GJ 581 d'           # Name of the planet
M_sun      = (0.31)*Msol          # Stellar mass. Inside( ) in solar masses. 
M_planet   = (7.1)*Mearth         # Mass of the planet. Inside ( ) in earth masses.
mu         = G*(M_sun+M_planet)
R          = (1.7)*Rearth         # Medium radius of the deformed body. Inside ( ) in earth radius. 
C          = 0.4*M_planet*R**2    # Inertia moment of the deformed body 
BmAC       = 5.0e-5               # Triaxiality [dimensionless]
rigidity   = 8.0e10               # Rigidity of the deformed body
tau        = 50*365.25*86400      # Andrade and Maxwell times 
alpha      = 0.2                  # Andrade's exponent [dimensionless]
e          = 0.27                 # Excentricity [dimensionless]
a          = (0.218)*AU           # Semi-major axis. Inside ( ) in astronomical units. 
#P          = 67.0*86400           # Orbital period
P          = np.sqrt(4*np.pi**2/mu * a**3)
#n          = 2.0*np.pi/P          # Mean motion
n          = np.sqrt(mu/a**3)
E0         = 0.0                  # Initial value for eccentric anomaly




############################################################
# Converting into canonical units

G        = 1
M_sun    = M_sun/uM       
M_planet = M_planet/uM
mu       = mu * uT**2/uL**3
R        = R/uL
C        = C/(uM*uL**2)
rigidity = rigidity*(uL*uT**2/uM)
tau      = tau/uT
a        = a/uL
P        = P/uT
n        = n*uT
E0       = 0.0


print "G        = ",G
print "M_sun    = ",M_sun
print "M_planet = ",M_planet
print "mu       = ",mu
print "R        = ",R
print "C        = ",C
print "rigidity = ",rigidity
print "tau      = ",tau
print "e        = ",e
print "a        = ",a
print "n        = ",n




############################################################
# A_l coefficients

def A_l(R,rigidity,M_planet,l):
    return 4*np.pi*R**4*rigidity*(2*l**2 + 4*l + 3)/(3*G*l*M_planet**2)





############################################################
# Tidal torque 

def tidal_torque(Omega,a,e):

    n = np.sqrt(mu/a**3)

    sum_tidal = 0.0
    for q in range(-1,8):
        omega = (2+q)*n - 2.0*Omega
        chi   = abs( omega )

        Re_compliance = 1 + ( (chi*tau)**(-alpha) )*np.cos(alpha*np.pi/2)*special.gamma(alpha+1.0)
        Im_compliance = -(chi*tau)**(-1.0) - ( (chi*tau)**(-alpha) )*np.sin(alpha*np.pi/2)*special.gamma(alpha+1.0)

        FR = -( 1.5*A_l(R,rigidity,M_planet,2)*Im_compliance ) / ( ( Re_compliance + A_l(R,rigidity,M_planet,2) )**2 + \
                                                                   Im_compliance**2 )*np.sign(omega)

        sum_tidal += ( G_function(e,2,0,q)**2 )*FR
    
    return 1.5*G*( M_sun**2 )*(R**5/a**6)*sum_tidal





############################################################
# Tidal torque times w_q

def tidal_torque_wq(Omega,a,e):
    
    n = np.sqrt(mu/a**3)

    sum_tidal_wq=0.0
    for q in range(-1,8):
        omega = (2+q)*n - 2.0*Omega
        chi   = abs( omega )

        Re_compliance = 1 + ( (chi*tau)**(-alpha) )*np.cos(alpha*np.pi/2)*special.gamma(alpha+1.0)
        Im_compliance = -(chi*tau)**(-1.0) - ( (chi*tau)**(-alpha) )*np.sin(alpha*np.pi/2)*special.gamma(alpha+1.0)

        FR = -( 1.5*A_l(R,rigidity,M_planet,2)*Im_compliance ) / ( ( Re_compliance + A_l(R,rigidity,M_planet,2) )**2 + \
                                                                   Im_compliance**2 )*np.sign(omega)

        sum_tidal_wq += ( G_function(e,2,0,q)**2 )*FR*omega

    return 1.5*G*( M_sun**2 )*(R**5/a**6)*sum_tidal_wq





############################################################
# Triaxial torque

def triaxial_torque(theta,a,e,t):

    n  = np.sqrt(mu/a**3)

    def kepler(E):
        return E-e*np.sin(E)-n*t
    E = newton(kepler,n*t)
    nu = 2.0*np.arctan( np.sqrt((1+e)/(1-e))*np.tan(0.5*E) )
    r=a*(1-e*np.cos(E))

    
    return -1.5*(C*BmAC)*(mu/a**3)*((a/r)**3)*np.sin(2.0*(theta-nu))    





############################################################
# Semimajor axis equation

def dadt(Omega,a,e):
    
    return -2*a**2/(G*M_sun*M_planet) * tidal_torque_wq(Omega,a,e) 





############################################################
# Eccentricity equation

def dedt(Omega,theta,a,e,t):
    
    return (tidal_torque(Omega,a,e)+triaxial_torque(theta,a,e,t)) * (a*(1-e**2)/(G*M_sun))**0.5 * 1.0/(a*e*M_planet) + dadt(Omega,a,e)/(2*a*e)*(1-e**2)





############################################################
# The angular velocity equation

def ang_velocity(Omega):
    return Omega




############################################################
# Integration parameters and initial conditions

max_dt = (100)*365.25*86400/uT # Maximum time step allowed. Inside ( ) in years
t_ini  = 0.0
t_end  = 50*P
#t_end = (1)*365.25*86400/uT
N      = 5000
time_array   = np.linspace(t_ini,t_end,N)


theta_ini = 0.0*(np.pi/180.0)
Omega_ini = 2.51*n 


initial_conditions = [theta_ini,Omega_ini,a,e] #Initial condition vector

print max_dt

# Printing info

print
print "Line number: ",N
print "t_end: %f yrs"%(t_end*uT/(86400*365.25))


def func(eta,t):
    theta = eta[0]
    Omega = eta[1]
    a     = eta[2]
    e     = eta[3]
    
    return [ang_velocity(Omega),
            (tidal_torque(Omega,a,e)+triaxial_torque(theta,a,e,t))/C,
            dadt(Omega,a,e),
            dedt(Omega,theta,a,e,t)]

print "Running ..."




############################################################
# The solution ...
start_time = time.time()
solucion,info = odeint(func,initial_conditions,time_array,full_output=True,printmessg=1)#hmax=max_dt)

print info['hu']


file1=open("evolution_corrected.dat","w")

for i in np.arange(0,len(time_array)):
    theta_int = solucion[i][0]
    Omega_int = solucion[i][1]
    a_int     = solucion[i][2]
    e_int     = solucion[i][3]
    n_int     = np.sqrt(mu/a_int**3)

   
    file1.write( "%1.5e     %1.5e    %1.5e    %1.9e    %1.9e   \n " % 
                 (time_array[i]*uT/(86400*365.25), theta_int, Omega_int/n_int, a_int, e_int) )
    
    """
    if (i%10==0.0):
        #(time   x   y   vx   vy   theta   omega/n   OrbitalEnergy)
        file1.write( "%1.5e     %1.5e    %1.5e     %1.5e    %1.5e   %1.5e   %1.5e   %1.5e \n " % 
                     (time[i]*uT,r[0],r[1],v[0],v[1],solucion[i][4],solucion[i][5]/n,OrbitalEnergy) )
    """


print
print "Execution time = %s hrs" %((time.time()-start_time)/3600.)
print
    


    
