"""In this module we compute the ray tracing in 2D based on
"Applied Acoustic" Qi Mo, Hengchin Yeh, Ming Lin, Dinesh Manocha  2016

WARNING : I am not sure that the precision of the roots solver is the same on
each computer reduce the precision if necessary in the function tracePeaceofRay.
If there is a problem some ray will pass through the ground but you can fix this
by putting a new starting when there is a reflexion a little bit more than me
over 0.

16/06/20:
Include reflexion on a ground f(x)=0
2D : plan of the axes of x and z

POSSIBLE IMPROVMENTS FOR A BETTER EFFIENCIENCY :
Resolve manually the intersections of the parabole and the circle or the ground.

Paralelle computing possible ?

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
import scipy.interpolate as inter


########################################
###HERE WE DEFINE FUNCTIONS USED AFTER###(see the body of the module just after)    
########################################

def traject(origin0,tetaOrigin0):
    """This function constructs the polynome which defines the trajectory
    of the peace of ray in an area where grad(V^-2) is considered to be constant.

    origin0 is the coordinates of the origin of the peace of ray in
    the coordinates system (x,y,z): (x0,y0,z0)
    
    tetaOrigin0 is the direction of the peace of ray at its origin : a scalar between -pi;pi"""

    #cf to "applied acoustic" to understand what follows in this function
    #the formulation of the ray traject in applied acoustic is true
    #for the local orthonormal coordinates (r,h) were h is the axis directed
    #by the direction of the gradient.
    
    x0,y0,z0=origin0
    alpha=norm(x0,y0,z0)
    Vorigin0=fV(x0,y0,z0) 
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
    """This function gives the ordinates corresponding to abscisses
    for a peace of ray from a origin within the range of validity
    and also return the new origin and direction.
    
    origin0 is the origin (first point) of the peace of ray
    direction0 is the direction at the origin0 in the (x,y,z) coordinate system"""

    x0,y0,z0=origin0
    xd,yd,zd=direction0
    validityRange=d(x0,y0,z0)
    tetaDirection=np.arctan(zd/xd)#we need an angle to compute the trajectory
    Z=traject(origin0,tetaDirection)
    
    x1=nextPoint(origin0,validityRange,Z)#x1 is the next by taking account of
    #the validity range only
    derivatZ=Z.deriv()#to get the new direction we will need to evaluate the
    #derivate of the trajectory at the next origin position


    ####first case if there is no possible reflexion on the ground####
    
    if z0-validityRange>0:#with this condition, it is impossible to encounter
    #the ground before de sphere of validity
    #so the next point is x1
        
        deriValue=derivatZ(x1)
        direction1=(1,0,deriValue)#the next direction vector
        
        z1=Z(x1)

        origin1=x1,y0,z1
        
    ####second case if there is a potential reflexion####
    
    else:#we can have a reflexion and we have to determine if it happens or not
        reflexion=False
        reflexionPoints=Z.r #it gives all the roots of the polynome describing
        #the peace of ray. In an other way it is all the intersection between 
        #the peace of ray and the ground
        
        for xsol in reflexionPoints:#maybe solutions are given in order and we dont need to browse the list
            if np.isreal(xsol) and xsol > x0 and xsol<=x1:#there will be reflexion
            #on the ground if there is a real intersection point which is
            #horizontally after the starting point of the peace of ray but
            #before x1 (abscisse of the intersection with the sphere of validity)
                x1=xsol.real#x1 is not the intersection with the circle anymore
                reflexion=True
        
                        
        if reflexion:
            
            deriValue=derivatZ(x1)
            direction1=(1,0,-deriValue)#the new vector of direction need to
            #symetric to the normal of the ground
            
            origin1=x1,y0,10**(-12) #we fix z to 0 because if we don't
            #do this, approximations can lead to a wrong reflexion
            #we are a little over 0 to evitate that approximation of the solver
            #leads to an other reflexion
            
            
        else:
            deriValue=derivatZ(x1)
            direction1=(1,0,deriValue)#the next direction vector

            z1=Z(x1)
            origin1=x1,y0,z1
            
    return(origin1,direction1,Z)
        

def nextPoint(origin0,validityRange,Z):
    """Return the next point in a free media. The next point is determine
    by the first intersection between the ray and the circle of validity"""
    x0,y0,z0=origin0
    
    #we want to resolve the abscisse of the intersection between a circle and a parabola
    #it leads to an 4th degrees equation whose coefficient are :

    a=Z[2]**2 #Z[i] is the coefficient of the i th degree of Z
    b=2*Z[2]*Z[1]
    c=Z[1]**2+2*Z[2]*(Z[0]-z0)+1
    d=2*Z[1]*(Z[0]-z0)-2*x0
    e=Z[0]**2-2*Z[0]*z0+z0**2+x0**2-validityRange**2
    
    validity_intersect=np.roots([a,b,c,d,e])
    
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


def gradOrientation(z0):
    """The axis z is normally oriented by the direction of grad(V^-2) but in
    this program we fix the z-axis as the vertical axe oriented upward. Assuming
    that the gradient direction is the z-axes (because in many cases media caracteritics
    depend on the high), we need to determine if the gradient is oriented upward
    or downward to finally correct the formula of the ray with a minus if necessary"""

    value=grad(0,0,z0)
    correction=np.sign(value)#correction =1 if value >0 and -1 if
    #value<0 (value<0 means that the gradient and the z axis are opposed)
    if correction==0:#Inthe case where the gradient is null we change nothing
        correction=1
    return (int(correction))


#########################
###Code to trace a ray###
#########################


#First we define ||grad(V^-2)|| as a numpy function which
#depends on the location. It will be used in other function

tck=inter.splrep(zProfil,cProfil,s=0)
fV=lambda x,y,z:np.float64(inter.splev(z,tck,ext=2))
V_2=cProfil**(-2)
tck2=inter.splrep(zProfil,V_2,s=0)
grad=lambda x,y,z:np.float64(inter.splev(z,tck2,ext=2,der=1))
norm=lambda x,y,z:np.abs(grad(x,y,z))
laplacian=lambda x,y,z:np.float64(inter.splev(z,tck,ext=2,der=2))

#norm is now a function of the space that returns
#||grad(V^-2(x,y,z))||for a given speed profile



#I changed the formula of delta given in applied acoustic :
#Now if the variation of the gradient is high relatively to its norm
#the range is low
#We had a second term in the case where the variation of the gradient is
#infinite in some points because without it the range could tend to 0
d=lambda x,y,z:delta*norm(x,y,z)/np.abs(laplacian(x,y,z))+delta/100*(1/(1+delta*norm(x,y,z)/np.abs(laplacian(x,y,z))))


def traceRay(rayon):
    
    origin0=rayon.origin
    direction0=rayon.direction
    i=0
    #we will construct  an array with four columns which for a ligne i gives
    #a,b,c,x1 (respectivly the first coefficient of the polynome i, the second,
    #the third and the last point of the peace of ray)
    #the first point can be deduced with the last "last point"
    
    #this form is not the more simple but enables to use the hdf5 saving

    
    ray=np.array([np.inf,np.inf,np.inf,origin0[0]])#for the first indice we just
    #need to save the origin point
    #coefficients can be used to save others thing like the size the ray for example 
    
    while origin0[0]<distance:####In this loop norm(x,y,z)
        #is computed 2 times (one time in the function d and one time in
        #the function traject) but I am not sure that it will improve performance
        #a lot. 
        origin1,direction1,Z=tracePeaceOfRay(origin0,direction0)
        
        ray=np.vstack((ray,np.hstack((np.array(Z),origin1[0]))))
        
        origin0=origin1 #we could have directly put origin0 instead in the first
        #line of the while 
        direction0=direction1

        i+=1 #to count the number of peace of ray
    return(ray,i)



    
    
    
