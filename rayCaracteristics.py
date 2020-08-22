""" This module enables to calculate some data about a ray """
import numpy as np
from Donnees import *

###Series of function that compute the caracteristic for a peace of ray###

def PeaceLength(a,b,c,x1,x0):
    """compute the length of a peace of ray"""
    X1=2*a*x1+b
    X0=2*a*x0+b
    peaceLength=((np.sqrt(X1**2+1)*X1+np.arcsinh(X1))-\
                (np.sqrt(X0**2+1)*X0+np.arcsinh(X0)))/(4*a)
    return(peaceLength)

def LocalApogee(Z,a,b,x1,x0):
    """find the apogee of a peace of ray"""
    Xa=-b/(2*a)#abscisse of the summit of the parabola
    if Xa>=x0 and Xa<x1:#if the summit is in the sphere of validity
        localApogee=Z(Xa)
    else:
        localApogee=np.maximum(Z(x1),Z(x0))#else the parabola is monotonic on the sphere
    return(localApogee)
    
     
def Peacett(Z,x1,x0,step,v):
    """compute travel time of a peace of ray"""
    peacett=0
    i1=x0
    for i2 in np.arange(x0+step,x1,step):
        t=np.sqrt(step**2+(Z(i2)-Z(i1))**2)/v(i2,Z(i2)) #the numerator is a ds (s the curvilign abscisse)
        peacett+=t
        i1=i2
    t=np.sqrt((x1-i1)**2+(Z(x1)-Z(i1))**2)/v(i1,Z(i1))#to integrate until x1
    peacett+=t
    return(peacett)

def LocalAttenuation(Z,x1,reflexion):
    if reflexion:
        derivatZ=Z.deriv()
        deriValue=derivatZ(x1)
        angle=np.arctan(deriValue)
        Zs=surfaceImpedance(impedance)#Zs is the normalised impedance of the ground
        LocalAttenuation=Rp(Zs,angle)
    else:
        LocalAttenuation=1
    return(LocalAttenuation)
def Rp(Zs,angle):
    return((Zs*np.cos(angle)-1)/(Zs*np.cos(angle)+1))#for plane waves
def surfaceImpedance(impedance):
    return(impedance)#locally reacting surface approximation

impedance=10 #possibility to implemente some code which compute
#a value of normalised impedance according to the surface caracteristic and frequency


###The following functions compute the caracteristics for a ray and then for
###the ray shooting


def datas(ray,v):
    """ To calculate the apogee, the ray length, the travel time of a ray and
the attenuation of a ray"""
    x0=ray[0,3]
    length=0
    step=distance/1000
    travelTime=0
    apogee=0
    length=0
    attenuation=1
    for a,b,c,x1,reflexion in ray[1:]:
        Z=np.poly1d([a,b,c])
        length+=PeaceLength(a,b,c,x1,x0)
        apogee=np.maximum(apogee,LocalApogee(Z,a,b,x1,x0))
        travelTime+=Peacett(Z,x1,x0,step,v)
        attenuation*=LocalAttenuation(Z,x1,reflexion)
        x0=x1
    return(length,apogee,travelTime,attenuation)

def RaysData(raysMatrix,raysIndex,v):
    raysData=np.empty(shape=(nbRay,4),dtype=np.float64)
    beginning=0
    for i,indice in enumerate(raysIndex):
        length,apogee,travelTime,attenuation=datas(raysMatrix[beginning:indice],v)
        raysData[i,]=[length,apogee,travelTime,attenuation]
        beginning=indice
    return(raysData)
        
    
