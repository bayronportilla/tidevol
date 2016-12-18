
############################################################
# Hansen coefficients calculator
# @bayronportilla-2016


from scipy import special
import numpy as np


############################################################
# Calcula coeficiente binomial con n<0 y k>=0

def binomial(n,k):
    return ((-1.0)**k)*special.binom(-n+k-1,k) # See for ref. https://arxiv.org/pdf/1105.3689.pdf


############################################################
# Calcula funcion G

def G_function(e,l,p,q):

    if(-2*l+2*p >= 0):
        s1 = -2*l+2*p
    elif(-2*l+2*p < 0):
        s1 = 100

    if(-2*p >= 0):
        t1 = -2*p
    elif(-2*p < 0):
        t1 = 100
        
    beta = e/(1+np.sqrt(1-e**2))
    
    sum = 0.0


    for s in range(0,s1+1):
        for t in range(0,t1+1):
            if(-2*l+2*p>=0 and -2*p>=0): 
                sum += special.binom(-2*l+2*p,s)*special.binom(-2*p,t)*(-beta)**(s+t) * special.jn(q-s+t,(l-2*p+q)*e)
            elif(-2*l+2*p>=0 and -2*p<0):
                sum += special.binom(-2*l+2*p,s)*binomial(-2*p,t)*(-beta)**(s+t) * special.jn(q-s+t,(l-2*p+q)*e)
            elif(-2*l+2*p<0 and -2*p>=0):
                sum += binomial(-2*l+2*p,s)*special.binom(-2*p,t)*(-beta)**(s+t) * special.jn(q-s+t,(l-2*p+q)*e)
            elif(-2*l+2*p<0 and -2*p<0):
                sum += binomial(-2*l+2*p,s)*binomial(-2*p,t)*(-beta)**(s+t) * special.jn(q-s+t,(l-2*p+q)*e)



    return (1 + beta**2)**l * sum





