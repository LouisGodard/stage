"""In this module we compute the ray tracing based on "Applied Acoustic"  """

from Donnees import *
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from sympy.vector import *
#sympy was already import in the datas module  


    

def traject(origin0,Vorigin0,tetaOrigin0):#origin0 to designe the orgin of
    #a peace of ray
    """This function built the function which defines the trajectory
of the ray in area where grad(V^-2) is considered to be constant.

origin has to have the form of (x,y,z)"""
    x0,y0,z0=origin0
    alpha=norm(x0,y0,z0)
    correction=int(gradOrientation(z0))
    teta=correction*tetaOrigin0#because if we change the orientation of z,
    #we also change angle direction
    
    epsilon=np.cos(teta)/Vorigin0
    zf=(Vorigin0**(-2)-epsilon**2)/alpha
    
    
    corrdirect=teta/abs(teta)#to correct a mistake in the
    #formula of rf (sign of teta matters for the summit position)
    rf=-corrdirect*2*epsilon*np.sqrt(Vorigin0**(-2)-epsilon**2)/alpha
    
    H=lambda r:(alpha*(r-rf)**2/(4*epsilon**2)-zf)*correction #H is the function of the
    #trajectory of a ray which pass through origin (0,0) with a direction teta
    #un moins a été mis dans un premier temps devant l'expression de la
    #trajectoire car z est orienté dans le même sens que grad(V-2)
    h=lambda r:H(r-x0)+z0 #h takes care that the origin is (x0,y0,z0) and not (0,0,0)
    return (h)


def tracePeaceOfRay(origin0,direction0):#origin and direction of the peace of ray
    """This function gives the ordinates corresponding to abscisses
for a peace of ray from a origin within the range of validity
and also return the new origin and direction.
origin is a triplet
direction is an angle"""
    x0,y0,z0=origin0
    xd,yd,zd=direction0
    validityRange=d(x0,y0,z0)
    tetaDirection=np.arctan(zd/xd)
    h=traject(origin0,fV(x0,y0,z0),tetaDirection)

    #In the following step we try to determine where to stop the tracing of the
    #peace of ray
    circle=(x-x0)**2+(z-z0)**2-validityRange**2
    #with this form circle=0 is the equation of a circle centered on (x0,z0)
    #with a validityRange radius length
    #we can find the intersections between our validity sphere and our ray by
    #resolving the system {z=h(x),circle=0} for x,z

    hInCircle=circle.subs(z,h(x)) #As the first step of resolving the system
    #we subsituted z by h(x) in circle
    #we could have used the function solve directly with the brut system but
    #by doing this we save computing time

    validity_intersect=solve(hInCircle,x)#here we determine
    # all the intersections of the ray and the validity shpere.
    #Be carefull the complex solutions are also returned


    #Now we will extract the point that will be the beginning of our new peace of ray
    x1=x0 #xnplusn1 will be the value of the next absciss origin
    for xsol in validity_intersect:
        if xsol.is_real and xsol>x0 and xsol>x1: #the next origin need to be real
            #upper than the actual origin
            #and we want to select the first/closest point that the ray intersects
            x1=xsol

    x1=float(x1)#convert to float because it is
    #important not to mix sympy and numpy
    #for exemple np.linspace does not support sympy.core.numbers.Float
    
    #to get the new direction we need to evaluate the derivate of the trajectory
    #at the new origin position
    derivative=diff(h(x),x) 
    deriValue=float(derivative.subs(x,x1))#convert to float because it is
    #important not to mix sympy and numpy
    direction1=(1,0,deriValue)#the vector of direction
    #Datas to trace the peace of ray
    absc=np.linspace(x0,x1,10)
    ordinates=h(absc)

    origin1=absc[-1],y,ordinates[-1]
    
    return(origin1,direction1,absc,ordinates)



def gradOrientation(z0):
    """The axis z is normally oriented by the direction of grad(V^-2) but in
this program we fix the z-axis as the vertical axe oriented upward. Assuming
that the gradient direction is the z-axes (because in many cases media caracteritics
depend on the high), we need to determine if the gradient is oriented upward
or downward to finally correct the formula of the ray with a minus if necessary"""

    projectionZ=R.k.projection(grad, scalar=True)
    value=projectionZ.subs(R.z,z0)
    correction=value/abs(value)#correction =1 if value >0 and -1 if
    #value<0 (which means that the gradient and the z axe are opposed
    return (correction)



#First we define ||grad(V^-2)|| as a numpy function which
#depends on the location. It will be used in other function
V_2=V**(-2) 
R = CoordSys3D('R')#R.z is the coordinate use by sympy for z
scaV_2=V_2.subs(z,R.z)
grad=gradient(scaV_2)#grad is used in gradOrientation
norm=grad.dot(grad)#grad.grad=||grad||^2
norm=norm**(1/2)
norm = lambdify((R.x,R.y,R.z), norm, "numpy")#this line convert norm to a lambda
#function from a symy expression

#norm is now a function of the space that returns
#||grad(V^-2(x,y,z))||for a given speed profile

delta=0.005
d=lambda x,y,z:(8*delta/norm(x,y,z))**(1/3)

abscisses=np.array([])
ordinates=np.array([])
distance=700
lastPointdistance=rayon.origin[0]
origin0=rayon.origin
direction0=rayon.direction
i=0
while lastPointdistance<distance and i<50:#####In this loop norm(x,y,z)
    #is computed 2 times (one time in the function d and one time in
    #the function traject) but I am not sure that it will improve performance
    #a lot. 
    origin1,direction1,absc,ordin=tracePeaceOfRay(origin0,direction0)
    abscisses=np.concatenate([abscisses,absc])
    ordinates=np.concatenate([ordinates,ordin])
    origin0=origin1 #we could have directly put origin0 instead in the first
    #line of the while 
    direction0=direction1
    lastPointdistance=origin1[0] #origin1[0] corresponds to x1
    i+=1

plt.plot(abscisses,ordinates)
sol=[0,distance]
plt.plot(sol,[0,0],"--")
plt.ylim(0,18)
plt.xlim(0,700)
plt.show()    



    
    
    
