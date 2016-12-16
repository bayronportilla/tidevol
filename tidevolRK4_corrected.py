#!/usr/bin/env python

from scipy.integrate import odeint,ode
from scipy.optimize import newton
from scipy.special import binom,gamma
from math import sqrt,cos,sin,tan,atan
from scipy import arctan,special
import numpy as np
from numpy.linalg import norm
from Units import units
import time
#import matplotlib.pyplot as plt
#import subprocess


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







#Calcula funcion de Bessel de orden v en la variable x
#print special.jn(v,x)


############################################################
# Hansen coefficients calculator

#Calcula binomiales con n<0 y k>=0
def binomial(n,k):
    return ((-1.0)**k)*special.binom(-n+k-1,k) # See for ex. https://arxiv.org/pdf/1105.3689.pdf


def G_function(e,q):
    beta = e/(1+sqrt(1-e**2))
    sum  = 0.0
    for s in range(0,100):
        sum += special.jn(q-s,(2+q)*e)*binomial(-4,s)*(-beta)**s
    return ((1 + beta**2)**2)*sum


def A2(R,rigidity,M_planet):
    return (38.0*np.pi*rigidity*R**4) / (3.0*G*M_planet**2) #Makarov's form [Makarov 2012]





############################################################
# Tidal torque only

def tidal_torque(Omega,a,e):

    n = np.sqrt(mu/a**3)
    sum_torque=0.0

    for q in range(-1,5):
        omega = (2+q)*n - 2.0*Omega
        chi   = abs( omega )

        Re_compliance = 1 + ( (chi*tau)**(-alpha) )*cos(alpha*np.pi/2)*special.gamma(alpha+1.0)
        Im_compliance = -(chi*tau)**(-1.0) - ( (chi*tau)**(-alpha) )*sin(alpha*np.pi/2)*special.gamma(alpha+1.0)

        FR = -( 1.5*A2(R,rigidity,M_planet)*Im_compliance ) / ( (Re_compliance + A2(R,rigidity,M_planet))**2 + Im_compliance**2 )*np.sign(omega)

        sum_torque += ( G_function(e,q)**2 )*FR
    
    return 1.5*G*( M_sun**2 )*(R**5/a**6)*sum_torque





############################################################
# Tidal torque times w_q

def tidal_torque_wq(Omega,a,e):
    
    sum_torque_wq=0.0
    n = sqrt(mu/a**3)

    for q in range(-1,5):
        omega = (2+q)*n - 2.0*Omega
        chi   = abs( omega )

        Re_compliance = 1 + ( (chi*tau)**(-alpha) )*cos(alpha*np.pi/2)*special.gamma(alpha+1.0)
        Im_compliance = -(chi*tau)**(-1.0) - ( (chi*tau)**(-alpha) )*sin(alpha*np.pi/2)*special.gamma(alpha+1.0)

        FR = -( 1.5*A2(R,rigidity,M_planet)*Im_compliance ) / ( (Re_compliance + A2(R,rigidity,M_planet))**2 + Im_compliance**2 )*np.sign(omega)

        sum_torque_wq += ( G_function(e,q)**2 )*FR*omega

    return 1.5*G*( M_sun**2 )*(R**5/a**6)*sum_torque_wq




############################################################
# Tidal+triaxial torque

def triaxial_torque(theta,a,e,t):

    n  = np.sqrt(mu/a**3)

    def kepler(E):
        return E-e*sin(E)-n*t
    E = newton(kepler,n*t)
    nu = 2.0*np.arctan( sqrt((1+e)/(1-e))*np.tan(0.5*E) )
    
    r=a*(1-e*np.cos(E))

    return -1.5*(C*BmAC)*(mu/a**3)*((a/r)**3)*sin(2.0*(theta-nu))    



############################################################
# Semimajor axis equation

def a_evol(Omega,a,e):
    return -2*a**2/(G*M_sun*M_planet) * tidal_torque_wq(Omega,a,e) 




############################################################
# Eccentricity equation

def e_evol(Omega,theta,a,e,t):
    return (tidal_torque(Omega,a,e)+triaxial_torque(theta,a,e,t)) * (a*(1-e**2)/(G*M_sun))**0.5 * 1.0/(a*e*M_planet) + a_evol(Omega,a,e)/(2*a*e)*(1-e**2)





############################################################
# The angular velocity equation

def ang_velocity(Omega):
    return Omega





############################################################
# Integration parameters and initial conditions

max_dt = (100)*365.25*86400/uT # Maximum time step allowed. Inside ( ) in years
t_ini  = 0.0
t_end  = 50000*P
#t_end = (1)*365.25*86400/uT
#N     = (t_end-t_ini)/h
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
            a_evol(Omega,a,e),
            e_evol(Omega,theta,a,e,t)]

print "Running ..."

start_time = time.time()


############################################################
# The solution ...

solucion,info = odeint(func,initial_conditions,time_array,full_output=True,printmessg=1)#hmax=max_dt)



print
print "execution time = %s hrs" %((time.time()-start_time)/3600.)
print

print info['hu']

exit()
"""
plt.figure()
plt.plot(time,solucion[:,1]/n)
plt.savefig("prueba.png")
"""

file1=open("evolution_corrected.dat","w")

for i in np.arange(0,len(time)):
    theta = solucion[i][0]
    Omega = solucion[i][1]
    """
    OrbitalEnergy=0.5*(norm(v))**2-mu/norm(r)
    """
    
    file1.write( "%1.5e     %1.5e    %1.5e   \n " % 
                 (time[i]*uT/(86400*365.25),theta,Omega/n) )
    
    """
    if (i%10==0.0):
        #(time   x   y   vx   vy   theta   omega/n   OrbitalEnergy)
        file1.write( "%1.5e     %1.5e    %1.5e     %1.5e    %1.5e   %1.5e   %1.5e   %1.5e \n " % 
                     (time[i]*uT,r[0],r[1],v[0],v[1],solucion[i][4],solucion[i][5]/n,OrbitalEnergy) )
    """

    


    
