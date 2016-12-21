#!/usr/bin/env python

############################################################
# Este modulo calcula el sistema de unidades canonicas
# apropiado para un problema determinado. Para usarlo,
# identifique primero dos unidades fundamentales que
# serviran para deducir la tercera. Supongamos que 
# fijamos uL y uM. Para determinar la unidad de tiempo
# simplemente haga units.(uL=xxxx, uM=xxxxx). El programa
# devolvera la unidad derivada en Sistema Internacional (SI)
# y las otras dos como un arreglo de la forma [uM,uL,uT],
# siempre en ese orden.

# ADVERTENCIA: el orden de los argumentos es irrelevante.
# ADVERTENCIA: las unidades ingresadas y las devueltas
# estan en SI. 

# @bayronportilla-2016

import numpy as np

############################################################
# Calcula sistema de unidades canonicas

G = 6.6740831e-11

def units(**kwargs):
    if( kwargs.keys()==['uM','uL'] or kwargs.keys()==['uL','uM'] ):
        uT = ( kwargs['uL']**3 / (G*kwargs['uM']) )**0.5
        return [kwargs.get('uM'), kwargs.get('uL'), uT]
    
    elif( kwargs.keys()==['uM','uT'] or kwargs.keys()==['uT','uM'] ):
        uL = ( G*kwargs['uM']*kwargs['uT']**2 )**(1.0/3.0)
        return [kwargs.get('uM'), uL, kwargs.get('uT')]
    
    elif( kwargs.keys()==['uT','uL'] or kwargs.keys()==['uL','uT'] ):
        uM = ( kwargs['uL']**3 / (kwargs['uT']**2 * G) )
        return [uM, kwargs.get('uL'), kwargs.get('uT')]






