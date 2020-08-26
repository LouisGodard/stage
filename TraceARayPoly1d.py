"""In this module we compute the ray tracing in 2D based on
"Analytic ray curve tracing for outdoor sound propagation"
Qi Mo, Hengchin Yeh, Ming Lin, Dinesh Manocha  2016, applied acoustic

WARNING : I am not sure that the precision of the roots solver is the same on
each computer reduce the precision if necessary in the function tracePeaceofRay.
If there is a problem some ray will pass through the ground but you can fix this
by putting a new starting when there is a reflexion a little bit more than me
over 0.

POSSIBLE IMPROVMENTS FOR A BETTER EFFIENCIENCY :
Use closed-form solution to resolve the intersection between the parabola and the circle.
Paralelle computing possible ?

x,z are the orthonormal cartesian coordinates where z is the vertical axis.
the suffixe 0 is assigned to a variable corresponding to the origin of a peace of ray.
the suffixe 1 is assigned to a variable corresponding to the last point
(origin of the next peace of ray) of a peace of ray.
direction is a direction given by a vector(xd,zd).
teta is a direction given by the angle between the horizontal axis and the direction vector
"""

from Donnees import *
import numpy as np
from sympy.vector import *
import scipy.interpolate as inter


#########################################
###HERE WE DEFINE FUNCTIONS USED AFTER###(see the body of the module just after)    
#########################################

def traject(origin0,tetaOrigin0):
    """This function constructs the polynome object which defines the trajectory
    of a peace of ray in an area where grad(V^-2) is considered constant.

    origin0 is the coordinates of the origin of the peace of ray in
    the coordinates system (x,z)
    
    tetaOrigin0 is the direction of the peace of ray at its origin : a scalar between -pi;pi"""

    #cf to "Analytic ray curve tracing for outdoor sound propagation"
    #to understand what follows in this function
    #the formulation of the ray traject in the article is true for
    #the local orthonormal coordinates (r,h) were h is the axis directed
    #by the direction of the gradient.
    
    x0,z0=origin0
    
    alpha=norm(x0,z0)
    Vorigin0=fV(x0,z0) 
    correction=gradOrientation(z0)#the local coordinates are oriented according
    #to the direction of the gradient of V^-2
    #if it is the opposite of the z axis then we will need to inverse the local
    #polynome
    teta=correction*tetaOrigin0#because changing the coordinates we inverted
    #the vertical axis orientation and so we also change angles orientation

    changeSign=np.sign(teta)#to correct a mistake in the
    #formula of rf (sign of teta matters for the vertex position)
    #In case where teta=0 it does not matter because rf is equal to 0


    #Now we calculated the local traject with the poly1d class
    #to do this I calculated the coefficient a,b,c of the polynom under this form:
    #ar^2+br+c
    #it appears that c=0 (the polynome pass throught the origin of the local coordinates)
    #I developped rf and epsilon to simplify the expression
    a=alpha/(4*(np.cos(teta)/Vorigin0)**2)
    b=changeSign*(1/np.cos(teta)**2-1)**(1/2)
    #In local coordinate it gives
    #h=correction*np.poly1d([a,b,0])
    #let's define a polynome for global coordinates
    #the polynome is Z(r)=h(r-x0)+z0 and by manual calculating we obtain:
    Z=correction*np.poly1d([a,b-2*a*x0,correction*z0-b*x0+a*x0**2])
    return (Z)


def tracePeaceOfRay(origin0,direction0):
    """This function obtains, with the function traject, the trajectory of the peace of ray whose origin
is origin0 and whose initial direction is direction 0
Then it calculates the limitating point due to the range
of validity or to an obstacle
It, finally, determines the initial direction for the next peace of ray
So it returns the traject Z,the limitating point : origin1
the next initial direction : direction1 and a 0 or a 1 depending on whether
there was a reflexion or not"""

    x0,z0=origin0
    xd,zd=direction0
    validityRange=d(x0,z0)
    tetaDirection=np.arctan(zd/xd)#we need an angle to compute the trajectory
    Z=traject(origin0,tetaDirection)
    x1=limitatingPoint(origin0,validityRange,Z)#x1 is the limitating point by taking
    #account of the validity range only
    derivatZ=Z.deriv()#to get the new direction we will need to evaluate the
    #derivate of the trajectory at the next origin position
    reflexion=0
    ####first case if there is no possible reflexion on the ground####
    
    if z0-validityRange>0:#with this condition, it is impossible to encounter
    #the ground before de sphere of validity
    #so the next point is x1
        
        deriValue=derivatZ(x1)
        direction1=(1,deriValue)#the next direction vector
        
        z1=Z(x1)

        origin1=x1,z1
        
    ####second case if there is a potential reflexion####
    
    else:#we can have a reflexion and we have to determine if it happens or not
        reflexionPoints=Z.r #it gives all the roots of the polynome describing
        #the peace of ray. In an other way it is all the intersections between 
        #the peace of ray and the ground
        
        for xsol in reflexionPoints:#maybe solutions are given in order and we dont need to browse the list
            if np.isreal(xsol) and xsol > x0 and xsol<=x1:#there will be reflexion
            #on the ground if there is a real intersection point which is
            #horizontally after the starting point of the peace of ray but
            #before x1 (abscisse of the intersection with the sphere of validity)
                x1=xsol.real#x1 is not the intersection with the circle anymore
                reflexion=1
                      
        if reflexion:
            
            deriValue=derivatZ(x1)
            direction1=(1,-deriValue)#the new vector of direction need to
            #symetric to the normal of the ground
            
            origin1=x1,10**(-10) #we fix z to 10**-10 because if we don't
            #do this, approximations can lead to a wrong reflexion
            #we are a little over 0 to evitate that approximation of the solver
            #leads to an other reflexion
            
            
        else:
            deriValue=derivatZ(x1)
            direction1=(1,deriValue)#the next direction vector

            z1=Z(x1)
            origin1=x1,z1
            
    return(origin1,direction1,Z,reflexion)
        

def limitatingPoint(origin0,validityRange,Z):
    """Return the next point in a free media. The next point is determine
    by the first intersection between the ray and the circle of validity"""
    x0,z0=origin0
    
    validity_intersect=solveIntersection(Z,origin0,validityRange)
    #we have determined
    #all the intersections of the ray and the validity shpere.
    #Be carefull the complex solutions are also returned

    #Now we will extract the point that will be the beginning of our new peace of ray
    x1=x0+validityRange #x1 can't be higher than x0 + validityRange
    for xsol in validity_intersect:
        if np.isreal(xsol) and xsol>=x0 and xsol<x1: #the next origin need to be real
            #upper than the actual origin
            #and we want to select the first/closest point that the ray intersects
            x1=xsol
    return(x1.real)

def solveIntersection(poly,center,radius):
    x0,z0=center
    #we want to resolve the abscisse of the intersection between a circle and a parabola
    #it leads to an 4th degrees equation whose coefficient are :

    a=poly[2]**2 #Z[i] is the coefficient of the i th degree of Z
    b=2*poly[2]*poly[1]
    c=poly[1]**2+2*poly[2]*(poly[0]-z0)+1
    d=2*poly[1]*(poly[0]-z0)-2*x0
    e=poly[0]**2-2*poly[0]*z0+z0**2+x0**2-radius**2
    
    return(np.roots([a,b,c,d,e]))
    

def gradOrientation(z0):
    """The axis z is normally oriented by the direction of grad(V^-2) but in
    this program we fix the z-axis as the vertical axe oriented upward. Assuming
    that the gradient direction is the z-axes (because in many cases media caracteritics
    depend on the high), we need to determine if the gradient is oriented upward
    or downward to finally correct the formula of the ray with a minus if necessary"""

    value=projectionZ(0,z0)
    correction=np.sign(value)#correction =1 if value >0 and -1 if
    #value<0 (value<0 means that the gradient and the z axis are opposed)
    if correction==0:#Inthe case where the gradient is null we change nothing
        correction=1
    return (int(correction))


#########################
###Code to trace a ray###
#########################


#First we treat the speed profile to create function that will be used by the rest of the program

if analytic:
    
    fV=lambdify((x,z),V,"numpy")#create a lambda function of V because it is more
    #efficient to evaluate the speed at a point with this
    
    V_2=V**(-2) 
    R = CoordSys3D('R')#R.z is the coordinate use by sympy for z
    scaV_2=V_2.subs(z,R.z)
    grad=gradient(scaV_2)#grad is used in gradOrientation

    projectionZ=R.k.projection(grad, scalar=True)#for the function gradOrientation
    projectionZ=lambdify((R.x,R.z), projectionZ, "numpy") #projectionZ is the projection of the grad on the z axis

    laplacian=divergence(grad)
    laplacian=lambdify((R.x,R.z), laplacian, "numpy")

    norm=grad.magnitude()
    norm = lambdify((R.x,R.z), norm, "numpy")#this line convert norm to a lambda
    #function from a symy expression

    #norm is now a function of the space that returns
    #||grad(V^-2(x,y,z))||for a given speed profile
else:
    tck=inter.splrep(zProfil,cProfil,s=0) #cubic interpolation
    #tck is 3-uple containing the knot-points,  , the coefficients  and the order
    #of the spline which 3 here
    fV=lambda x,z:np.float64(inter.splev(z,tck,ext=2)) #splev enables to evaluate the interpolation at a point
    V_2=cProfil**(-2)
    tck2=inter.splrep(zProfil,V_2,s=0) #interpolation of V^-2
    projectionZ=lambda x,z:np.float64(inter.splev(z,tck2,ext=2,der=1))#projectionZ is the projection of the grad on the z axis
    #this is the derivative but it is also the values of the grad because
    #we suppose that speed only variate with the altitude
    norm=lambda x,z:np.abs(projectionZ(x,z))
    laplacian=lambda x,z:np.float64(inter.splev(z,tck2,ext=2,der=2))

#I changed the formula of delta given in applied acoustic :
#Now if the variation of the gradient is high relatively to its norm
#the range is low
#We had a second term in the case where the variation of the gradient is
#infinite in some points because without it the range could tend to 0
#d=lambda x,z:delta*norm(x,z)/np.abs(laplacian(x,z))+delta*(1/(1+delta*norm(x,z)/np.abs(laplacian(x,z))))
d=lambda x,z:(8*delta/norm(x,z))**(1/3)


def traceRay(rayon):
    
    origin0=rayon.origin
    direction0=rayon.direction
    i=0
    #we will construct  an array with four columns which for a ligne i gives
    #a,b,c,x1 (respectivly the first coefficient of the polynome i, the second,
    #the third and the last point of the peace of ray)
    #the first point can be deduced with the last "last point"
    
    #this form is not the more simple but enables to use the hdf5 saving

    
    ray=np.array([np.nan,np.nan,np.nan,origin0[0],0])#for the first indice we just
    #need to save the origin point
    #coefficients can be used to save others thing like the size the ray for example 
    
    while origin0[0]<distance:####In this loop norm(x,y,z)
        #is computed 2 times (one time in the function d and one time in
        #the function traject) but I am not sure that it will improve performance
        #a lot. 
        origin1,direction1,Z,reflexion=tracePeaceOfRay(origin0,direction0)
        
        ray=np.vstack((ray,np.hstack((np.array(Z),origin1[0],reflexion))))
        
        origin0=origin1 #we could have directly put origin0 instead in the first
        #line of the while 
        direction0=direction1

        i+=1 #to count the number of peace of ray
    
    return(ray,i)



    
    
    
