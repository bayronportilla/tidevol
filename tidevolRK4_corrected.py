#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import newton
from scipy.special import binom,gamma
from scipy import arctan,special
from numpy.linalg import norm
from Units import units
from Ecc_Polynomial import G_function
import time
import params


############################################################
# Constants and definitions. All of it in SI

G       = 6.67408e-11
M_sun   = 1.988500e30 
M_earth = 5.9723e24 
R_earth = 6371000.0 
AU      = 149.6e9
InRad   = np.pi/180
InDeg   = 180/np.pi 


############################################################
# Canonical units

Unidades_Canonicas = units(uM=0.001*M_earth,uL=0.218*AU)
uM = Unidades_Canonicas[0]
uL = Unidades_Canonicas[1]
uT = Unidades_Canonicas[2]



############################################################
# Loading data from params file and converting to SI.

name       = params.name
M_s        = params.M_s * M_sun
M_p        = params.M_p * M_earth
mu         = G*(M_s+M_p) 
R          = params.R * R_earth
C          = 0.4*M_p*R**2 
BmAC       = params.BmAC
rigidity   = params.rigidity
tau        = params.tau * 365.25*86400
alpha      = params.alpha
e          = params.e
a          = params.a * AU
#P          = 67.0*86400           # Orbital period
P          = np.sqrt(4*np.pi**2/mu * a**3)
#n          = 2.0*np.pi/P          # Mean motion
n          = np.sqrt(mu/a**3)
E0         = params.E0 * InRad

t_ini      = params.t_ini * 365.25*86400
t_end      = params.t_end * 365.25*86400

theta_ini  = params.theta_ini * InRad
Omega_ini  = params.p * n 
a_ini      = params.a_ini * AU
e_ini      = params.e_ini




############################################################
# Converting into canonical units

G         = 1
M_s       = M_s/uM       
M_p       = M_p/uM
mu        = mu * uT**2/uL**3
R         = R/uL
C         = C/(uM*uL**2)
rigidity  = rigidity*(uL*uT**2/uM)
tau       = tau/uT
a         = a/uL
P         = P/uT
n         = n*uT
E0        = 0.0

t_ini     = t_ini/uT
t_end     = t_end/uT

theta_ini = theta_ini
Omega_ini = Omega_ini * uT
a_ini     = a_ini/uL
e_ini     = e_ini


print "G        = ",G
print "M_s      = ",M_s
print "M_p      = ",M_p
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

def A_l(R,rigidity,M_p,l):
    return 4*np.pi*R**4*rigidity*(2*l**2 + 4*l + 3)/(3*G*l*M_p**2)





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

        FR = -( 1.5*A_l(R,rigidity,M_p,2)*Im_compliance ) / ( ( Re_compliance + A_l(R,rigidity,M_p,2) )**2 + \
                                                                   Im_compliance**2 )*np.sign(omega)

        sum_tidal += ( G_function(e,2,0,q)**2 )*FR
    
    return 1.5*G*( M_s**2 )*(R**5/a**6)*sum_tidal





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

        FR = -( 1.5*A_l(R,rigidity,M_p,2)*Im_compliance ) / ( ( Re_compliance + A_l(R,rigidity,M_p,2) )**2 + \
                                                                   Im_compliance**2 )*np.sign(omega)

        sum_tidal_wq += ( G_function(e,2,0,q)**2 )*FR*omega

    return 1.5*G*( M_s**2 )*(R**5/a**6)*sum_tidal_wq





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
# The angular velocity equation

def dthetadt(Omega):
    return Omega





############################################################
# Semimajor axis equation

def dadt(Omega,a,e):
    
    return -2*a**2/(G*M_s*M_p) * tidal_torque_wq(Omega,a,e) 




############################################################
# Eccentricity equation

def dedt(Omega,theta,a,e,t):
    
    return (tidal_torque(Omega,a,e)+triaxial_torque(theta,a,e,t)) * (a*(1-e**2)/(G*M_s))**0.5 * 1.0/(a*e*M_p) + dadt(Omega,a,e)/(2*a*e)*(1-e**2)






############################################################
# Initial conditions, time array, and output line number 

eta_ini      = [theta_ini,Omega_ini,a_ini,e_ini] 
N            = params.N
time_array   = np.linspace(t_ini,t_end,N)

max_dt = (100)*365.25*86400/uT # Maximum time step allowed. Inside ( ) in years

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
    
    return [dthetadt(Omega),
            (tidal_torque(Omega,a,e)+triaxial_torque(theta,a,e,t))/C,
            dadt(Omega,a,e),
            dedt(Omega,theta,a,e,t)]

print "Running ..."




############################################################
# The solution ...

start_time = time.time()
solucion,info = odeint(func,eta_ini,time_array,full_output=True,printmessg=1)

print info['hu']
exit()



############################################################
# Saving info ...

file1=open("evolution_corrected.dat","w")

for i in np.arange(0,len(time_array)):
    theta_int = solucion[i][0]
    Omega_int = solucion[i][1]
    a_int     = solucion[i][2]
    e_int     = solucion[i][3]
    n_int     = np.sqrt(mu/a_int**3)

   
    file1.write( "%1.5e     %1.5e    %1.5e    %1.9e    %1.9e   \n " % 
                 (time_array[i]*uT/(86400*365.25), theta_int, Omega_int/n_int, a_int, e_int) )
 

print
print "Execution time = %s hrs" %((time.time()-start_time)/3600.)
print
    


    
