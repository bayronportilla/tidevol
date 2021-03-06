#!/usr/bin/env python

import numpy as np
from scipy.integrate import odeint,ode
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
tau_M      = params.tau_M * 365.25*86400
alpha      = params.alpha
e          = params.e
a          = params.a * AU
P          = params.P * 86400         
n          = 2.0*np.pi/P 
E0         = params.E0 * InRad

t_ini      = params.t_ini * 365.25*86400
t_end      = params.t_end * 365.25*86400

theta_ini  = params.theta_ini * InRad
Omega_ini  = params.p * n 
a_ini      = params.a_ini * AU
e_ini      = params.e_ini



############################################################
# Printing info ...

print
print "Bulk properties:"
print "----------------"
print "Object             :",name
print "Stellar mass       = %.3f        [Solar masses]"  % (M_s/M_sun)
print "Planetary mass     = %.3f        [Earth masses]"  % (M_p/M_earth)
print "Planetary radius   = %.3f        [Earth radius]"  % (R/R_earth)
print "Triaxiality        = %.3e    [Dimensionless]"     % (BmAC)
print "Unrelaxed rigidity = %.3e    [Pascals]"           % (rigidity)
print "Relaxation time    = %.3f       [years]"          % (tau_M/(365.25*86400))
print "Andrade exponent   = %.3f        [Dimensionless]" % (alpha)

print
print "Dynamical properties:"
print "---------------------"
print "Eccentricity       = %.3f        [Dimensionless]" % (e)
print "Semimajor axis     = %.3f        [AU]"            % (a/AU)
print "Orbital period     = %.3f       [Days]"          % (P/86400)

print
print "Initial conditions and integration parameters:"
print "----------------------------------------------"
print "Initial sidereal angle    = %.3f      [Degrees]"  % (theta_ini*InDeg)
print "Initial rotational period = %.3f     [Days]"      % (2*np.pi/(params.p * n)/86400)
print "t_ini                     = %.3f      [years]"    % (t_ini/(365.26*86400))
print "t_end                     = %.3f      [years]"    % (t_end/(365.26*86400))
print "Number of output data     = %d               "    % params.N
print 





############################################################
# Writting input information in info.log file ...

file_3=open("info_%s.log"%(name),"w")
file_3.write("""
Bulk properties: 
---------------- 
Object             : %s
Stellar mass       = %.3f        [Solar masses]
Planetary mass     = %.3f        [Earth masses]  
Planetary radius   = %.3f        [Earth radius]  
Triaxiality        = %.3e        [Dimensionless]   
Unrelaxed rigidity = %.3e        [Pascals]       
Relaxation time    = %.3f        [years]          
Andrade exponent   = %.3f        [Dimensionless]


Dynamical properties:
---------------------
Eccentricity       = %.3f        [Dimensionless]
Semimajor axis     = %.3f        [AU]            
Orbital period     = %.3f        [Days]          


Initial conditions and integration parameters:
----------------------------------------------
Initial sidereal angle    = %.3f      [Degrees] 
Initial rotational period = %.3f      [Days]     
t_ini                     = %.3f      [years]   
t_end                     = %.3f      [years]   
Number of output data     = %d                
"""%(name,(M_s/M_sun),(M_p/M_earth),(R/R_earth),(BmAC),(rigidity),(tau_M/(365.25*86400)),(alpha),
     (e),(a/AU),(P/86400),(theta_ini*InDeg),(2*np.pi/(params.p * n)/86400), (t_ini/(365.26*86400)),
     (t_end/(365.26*86400)), params.N))






############################################################
# Converting into canonical units

G         = 1
M_s       = M_s/uM       
M_p       = M_p/uM
mu        = mu * uT**2/uL**3
R         = R/uL
C         = C/(uM*uL**2)
rigidity  = rigidity*(uL*uT**2/uM)
tau_M     = tau_M/uT
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





############################################################
# A_l coefficients

def A_l(R,rigidity,M_p,l):
    return 4*np.pi*R**4*rigidity*(2*l**2 + 4*l + 3)/(3*G*l*M_p**2)





############################################################
# Tidal torque 

def tidal_torque(Omega,a,e):

    sum_tidal = 0.0
    for q in range(-1,8):
        omega = (2+q)*n - 2.0*Omega
        chi   = abs( omega )
        tau_A = (100.0*np.exp(-chi/0.2) + 1.0)*tau_M

        Re_compliance = 1 + ( (chi*tau_A)**(-alpha) )*np.cos(alpha*np.pi/2)*special.gamma(alpha+1.0)
        Im_compliance = -(chi*tau_M)**(-1.0) - ( (chi*tau_A)**(-alpha) )*np.sin(alpha*np.pi/2)*special.gamma(alpha+1.0)

        FR = -( 1.5*A_l(R,rigidity,M_p,2)*Im_compliance ) / ( ( Re_compliance + A_l(R,rigidity,M_p,2) )**2 + \
                                                                   Im_compliance**2 )*np.sign(omega)

        sum_tidal += ( G_function(e,2,0,q)**2 )*FR
    
    return 1.5*G*( M_s**2 )*(R**5/a**6)*sum_tidal





############################################################
# Tidal torque times w_q

def tidal_torque_wq(Omega,a,e):
    
    sum_tidal_wq=0.0
    for q in range(-1,8):
        omega = (2+q)*n - 2.0*Omega
        chi   = abs( omega )
        tau_A = (100.0*np.exp(-chi/0.2) + 1.0)*tau_M

        Re_compliance = 1 + ( (chi*tau_A)**(-alpha) )*np.cos(alpha*np.pi/2)*special.gamma(alpha+1.0)
        Im_compliance = -(chi*tau_M)**(-1.0) - ( (chi*tau_A)**(-alpha) )*np.sin(alpha*np.pi/2)*special.gamma(alpha+1.0)

        FR = -( 1.5*A_l(R,rigidity,M_p,2)*Im_compliance ) / ( ( Re_compliance + A_l(R,rigidity,M_p,2) )**2 + \
                                                                   Im_compliance**2 )*np.sign(omega)

        sum_tidal_wq += ( G_function(e,2,0,q)**2 )*FR*omega

    return 1.5*G*( M_s**2 )*(R**5/a**6)*sum_tidal_wq




############################################################
# Triaxial torque

def triaxial_torque(theta,a,e,t):

    
    def kepler(E):
        return E-e*np.sin(E)-n*t

    E = newton(kepler,n*t)

    nu = 2.0*np.arctan2(1 , 1/(np.sqrt((1+e)/(1-e))*np.tan(0.5*E)) )
    
    r = a*(1-e*np.cos(E))
    
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

N            = params.N
eta_ini      = [theta_ini,Omega_ini] 
time_array   = np.linspace(t_ini,t_end,N)


max_dt = (0.0015)*365.25*86400/uT # Maximum time step allowed. Inside ( ) in years


def func(t,eta):
    theta = eta[0]
    Omega = eta[1]
    
    return [dthetadt(Omega),
            (tidal_torque(Omega,a,e)+triaxial_torque(theta,a,e,t))/C]

#print "Status: running ..."





############################################################
# The solution ...

start_time = time.time()
solucion = np.zeros([int(N),3])

solver = ode(func).set_integrator("dop853")#,method="bdf")
solver.set_initial_value(eta_ini)
for i in range(0,int(N)):
    solucion[i][0] = time_array[i]
    solucion[i][1] = solver.integrate(time_array[i])[0]
    solucion[i][2] = solver.integrate(time_array[i])[1]
    


#solucion,info = odeint(func,eta_ini,time_array,full_output=True,printmessg=1,hmax=max_dt,mxstep=1000)

#print info['hu']






############################################################
# Saving info ...

file_1 = open("evolution_fast_%s.dat"%(name),"w")
#file_2 = open("step_size_fast.log","w")

for i in range(0,len(time_array)):
    theta_int = solucion[i][1]
    Omega_int = solucion[i][2]
       
    file_1.write( "%1.5e     %1.5e    %1.5e \n " % 
                 (time_array[i]*uT/(86400*365.25), theta_int, Omega_int/n) )
 
"""
for i in range(0,len(info['hu'])):
    file_2.write("%d   %.17e \n" % (i,info['hu'][i]*uT/86400))
"""    

print
print "Execution time = %s hrs" %((time.time()-start_time)/3600.)
print
    


    
