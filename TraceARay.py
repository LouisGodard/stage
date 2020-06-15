"""In this module we compute the ray tracing in 2D based on
"Applied Acoustic" Qi Mo, Hengchin Yeh, Ming Lin, Dinesh Manocha  2016

WARNING : I am not sure that the precision of the sympy solver is the same on
each computer reduce the precision if necessary in the function tracePeaceofRay

15/06/20:
Include reflexion on a ground f(x)=0
2D : plan of the axes of x and z

x,y,z are the orthonormal cartesian coordinates where z is the vertical axis.
the suffixe 0 is assigned to a variable corresponding to the first point
(origin) of a pay of ray.
the suffixe 1 is assigned to a variable corresponding to the last point
(origin of the next peace of ray) of a peace of ray.
direction is a direction given by a vector(xd,yd,zd).
teta is a direction given by the angle between the horizontal axis and the vector of direction
origin is the cartesian coordinates of the origins points."""

from Donnees import *
import numpy as np
#import scipy as sp
import matplotlib.pyplot as plt
from sympy.vector import *
#sympy was already import in the datas module  


    

def traject(origin0,Vorigin0,tetaOrigin0):
    """This function built the function which defines the trajectory
    of the peace of ray in an area where grad(V^-2) is considered to be constant.

    origin0 is the coordinates of the origin of the peace of ray in
    the coordinates system (x,y,z): (x0,y0,z0)
    Vorigin0 is the sound speed at the origin : a scalar
    tetaOrigin0 is the direction of the peace of ray at its origin : a scalar between -pi;pi"""

    #cf to "applied acoustic" to understand what follows in this function
    #the formulation of the ray traject in applied acoustic is true
    #for the local orthonormal coordinates (r,h) were h is the axis directed
    #by the direction of the gradient.
    
    x0,y0,z0=origin0
    alpha=norm(x0,y0,z0)
    correction=gradOrientation(z0)
    teta=correction*tetaOrigin0#because changing the coordinates we inverted
    #the vertical axis orientation and so we also change angles orientation
    
    epsilon=np.cos(teta)/Vorigin0
    zf=(Vorigin0**(-2)-epsilon**2)/alpha
    
    
    changeSign=np.sign(teta)#to correct a mistake in the
    #formula of rf (sign of teta matters for the vertex position)
    #In case where teta=0 it does not matter because rf is equal to 0
    rf=-changeSign*2*epsilon*np.sqrt(Vorigin0**(-2)-epsilon**2)/alpha

    h=lambda r:(alpha*(r-rf)**2/(4*epsilon**2)-zf)#in local coordinates
    
    z=lambda r:correction*h(r-x0)+z0 #we change the coordinates system taking
    #to account that the vertical axis is inverted and that the origin point of 
    #the peace of ray is now (x0,y0,z0) in the (x,y,z) coordinates
    return (z)


def tracePeaceOfRay(origin0,direction0):#origin and direction of the peace of ray
    """This function gives the ordinates corresponding to abscisses
    for a peace of ray from a origin within the range of validity
    and also return the new origin and direction.
    
    origin0 is the origin (first point) of the peace of ray
    direction0 is the direction at the origin0 in the (x,y,z) coordinate system"""

    x0,y0,z0=origin0
    xd,yd,zd=direction0
    validityRange=d(x0,y0,z0)
    tetaDirection=np.arctan(zd/xd)
    Z=traject(origin0,fV(x0,y0,z0),tetaDirection)
    #we have the trajectory of the peace of ray, now we need to determine
    #where to stop it

    x1=nextPoint(origin0,validityRange,Z)#here x1 is finding the constraint
    #of the validity range only
    
    derivative=diff(Z(x),x)#to get the new direction we will need to evaluate the
    #derivate of the trajectory at the next origin position

    
    ####first case if there is no possible reflexion on the ground####
    
    if z0-validityRange>0:#with this condition, it is impossible to encounter
    #the ground before de sphere of validity
        deriValue=float(derivative.subs(x,x1))#convert to float because it is
        #important not to mix sympy and numpy
        direction1=(1,0,deriValue)#the next direction vector

        x1=float(x1)#convert to float because it is
        #important not to mix sympy and numpy
        #for exemple np.linspace does not support sympy.core.numbers.Float
    
        #Datas to trace the peace of ray
        absc=np.linspace(x0,x1,nbPoint)
        ordinates=Z(absc)

        origin1=absc[-1],y,ordinates[-1]
        
    ####second case if there is a potential reflexion####
    
    else:#we can have a reflexion and we have to determine if it happens or not
        reflexion=False
        reflexionPoints=solve(Z(x),x) #we determine all the intersection between
        #the peace of ray and the ground. The solver have a precision
        #of 10^-13 which can lead to errors that's why we add a 10**(-11) under
        
        for xsol in reflexionPoints:#maybe solution are given in order and we dont need to browse the list
            if xsol.is_real and xsol >= x0+10**(-11) and xsol<=x1:#there will be reflexion
            #on the ground if there is a real intersection point which is
            #horizontally after the starting point of the peace of ray but
            #before x1 (abcissa the intersection with the sphere of validity)
                x1=xsol#x1 is not the intersection with the circle anymore
                reflexion=True
                
        if reflexion:
            deriValue=float(derivative.subs(x,x1))#convert to float because it is
            #important not to mix sympy and numpy
            direction1=(1,0,-deriValue)#the new vector of direction need to
            #symetric to the normal of the ground

            x1=float(x1)#convert to float because it is
            #important not to mix sympy and numpy
            #for exemple np.linspace does not support sympy.core.numbers.Float
            
            #Datas to trace the peace of ray
            absc=np.linspace(x0,x1,nbPoint)
            ordinates=Z(absc)
            origin1=absc[-1],y,ordinates[-1] #we fix z to 0 because if we don't do this
            #approximations can lead to a wrong reflexion
            
        else:
            deriValue=float(derivative.subs(x,x1))#convert to float because it is
            #important not to mix sympy and numpy
            direction1=(1,0,deriValue)#the next direction vector

            x1=float(x1)#convert to float because it is
            #important not to mix sympy and numpy
            #for exemple np.linspace does not support sympy.core.numbers.Float
    
            #Datas to trace the peace of ray
            absc=np.linspace(x0,x1,nbPoint)
            ordinates=Z(absc)
            origin1=absc[-1],y,ordinates[-1]
    
    return(origin1,direction1,absc,ordinates)
        

def nextPoint(origin0,validityRange,Z):
    """Return the next point in a free media. The next point is determine
    by the first intersection between the ray and the circle of validity"""
    x0,y0,z0=origin0
    #In the following step we try to determine where to stop the tracing of the
    #peace of ray
    circle=(x-x0)**2+(z-z0)**2-validityRange**2
    #with this form circle=0 is the equation of a circle centered on (x0,z0)
    #with a validityRange radius length
    #we can find the intersections between our validity sphere and our ray by
    #resolving the system {z=h(x),circle=0} for x,z

    ZInCircle=circle.subs(z,Z(x)) #As the first step of resolving the system
    #we subsituted z by h(x) in circle
    #we could have used the function solve directly with the brut system but
    #by doing this we save computing time

    validity_intersect=solve(ZInCircle,x)#here we determine
    # all the intersections of the ray and the validity shpere.
    #Be carefull the complex solutions are also returned


    #Now we will extract the point that will be the beginning of our new peace of ray
    x1=x0 #xnplusn1 will be the value of the next absciss origin
    for xsol in validity_intersect:#maybe solution are given in order and we dont need to browse the list
        if xsol.is_real and xsol>x0 and xsol>x1: #the next origin need to be real
            #upper than the actual origin
            #and we want to select the first/closest point that the ray intersects
            x1=xsol
    return(x1)


def gradOrientation(z0):
    """The axis z is normally oriented by the direction of grad(V^-2) but in
    this program we fix the z-axis as the vertical axe oriented upward. Assuming
    that the gradient direction is the z-axes (because in many cases media caracteritics
    depend on the high), we need to determine if the gradient is oriented upward
    or downward to finally correct the formula of the ray with a minus if necessary"""

    projectionZ=R.k.projection(grad, scalar=True)
    value=projectionZ.subs(R.z,z0)
    correction=value/abs(value)#correction =1 if value >0 and -1 if
    #value<0 (value<0 means that the gradient and the z axis are opposed)
    correction=int(correction)#for sympy/numpy comptibility
    return (correction)



#First we define ||grad(V^-2)|| as a numpy function which
#depends on the location. It will be used in other function
V_2=V**(-2) 
R = CoordSys3D('R')#R.z is the coordinate use by sympy for z
scaV_2=V_2.subs(z,R.z)
grad=gradient(scaV_2)#grad is used in gradOrientation
norm=grad.magnitude()
norm = lambdify((R.x,R.y,R.z), norm, "numpy")#this line convert norm to a lambda
#function from a symy expression

#norm is now a function of the space that returns
#||grad(V^-2(x,y,z))||for a given speed profile

delta=0.01
d=lambda x,y,z:(8*delta/norm(x,y,z))**(1/3)
distance=700
nbPeace=distance/d(0,0,0)#gives approxmatly the number of peace of ray
nbPoint=int(50/nbPeace)#number of point to trace a peace of ray to have approximatly
#50 points to trace the trajectory of the ray

def traceRay():
    abscisses=np.array([])
    ordinates=np.array([])
    
    lastPointdistance=rayon.origin[0]
    origin0=rayon.origin
    direction0=rayon.direction
    i=0
    while lastPointdistance<distance and i<50:#####In this loop norm(x,y,z)
        #is computed 2 times (one time in the function d and one time in
        #the function traject) but I am not sure that it will improve performance
        #a lot. 
        NextOrigin,NextDirection,absc,ordin=tracePeaceOfRay(origin0,direction0)
        abscisses=np.concatenate([abscisses,absc])
        ordinates=np.concatenate([ordinates,ordin])
        origin0=NextOrigin #we could have directly put origin0 instead in the first
        #line of the while 
        direction0=NextDirection
        lastPointdistance=NextOrigin[0] #NextOrigin[0] corresponds to x1
        i+=1

    plt.plot(abscisses,ordinates)
    sol=[0,distance]
    plt.plot(sol,[0,0],"--")
    plt.ylim(-2,18)
    plt.xlim(0,700)
    plt.show()

import time
def no(x):
    t0=time.time()
    for i in range(1000):
        s=np.sign(x)
        if s==0:
            s=1
    t1=time.time()
    return(t1-t0)
def no2(x):
    t0=time.time()
    for i in range(1000):
        s=x/np.abs(x)
    t1=time.time()
    return(t1-t0)

    
    
    
