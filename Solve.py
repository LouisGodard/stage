"""In this module we develop mathematic function especially for
our algorithm"""

def intersectCircleParabol(radius, center, a, b, c):
    """Radius and center define the circle and a,b,c is are the coefficients of the
parabola. Center is only a scalar x0 because in our programm the center is on
the parabola (so z0=f(x))
To return the intersections, we will use the common technics to solve
quartic equation"""

    #making the changment of variable X=x-x0 and Z=z-z0
    #the resolution of the system lead to resolve
    #a2X4+2a(2x0a+b)X3+(1+(2x0a+b)2)X2-raidus2=0
    #this step was not necessary but make things easier
    #we set:
    #A=a**2
    #B=2*a*(2*x0*a+b)
    #C=1+(2*x0*a+b)**2
    #E=-radius

    #now we eliminate the third degree coefficient with
    #the changment of variable X=y-b/(4a) and we divided by a
    #It gives:
    #y4+py2+qy+r=0
    #with: 
    #p=C/A-3*B**2/(8*A**2)
    #q=-BC/(2*A**2)+B**3/(8*a**3)
    #r=E/A+CB**2/(16*A**3)-3*B**4/(256*A**4)

    #we need now to resolve
    #z**3+2*p*z**2+(p^2-4*r)z-q**2

    #we need to eliminate the second degree term with the changement
    #of variable (z=Y-3/(3*1)
    #It gives:
    #Y3+PY+Q=0
    #with:
    #P=(p^2-4*r)-(2*p)**2/3
    #Q=q**2-(p^2-4*r)*2*p/3+2/27*(2*p)**3

    




def coefficient(correction,alpha,epsilon,rf,zf,x0,validityRange):
    """WARNING THIS FUNCTION DOESNT WORK
This function return the coefficient of the quartic equation
that enables to find the intersection point between the parabola and the circle
to use the roots function of numpy
the system is :
(x-x0)^2+(z-z0)^2=validityRange
z=z(x)

the z0 in z(x) eliminate the one in the first equation so we dont need z0
in the arguments"""
    #by calcul
    a=correction**2*alpha**2/(16*epsilon**4)
    b=-4*(x0+rf)*alpha/(16*epsilon**4)*correction**2
    c=(-alpha*zf/(2*epsilon**2)+\
       6*(x0+rf)**2*alpha**2/(16*epsilon**4))*correction**2\
       +1
    d=((x0+rf)*alpha/epsilon**2*zf-4*(x0+rf)**3*alpha**2/(16*epsilon**4))*correction**2\
       -2*x0
    e=correction**2*(\
        zf**2-alpha*zf/(2*epsilon**2)*(x0+rf)**2+alpha**2/(16*epsilon**4)*(x0+rf)**4)\
        +x0**2-validityRange**2
    return(a,b,c,d,e)

def coefficient2(correction,alpha,epsilon,rf,zf,validityRange):
    """ This function return the coefficient of the quartic equation
that enables to find the intersection point between the parabola and the circle
in THE LOCAL COORDINATE SYSTEM
to use the roots function of numpy
the system is :
(x-x0)^2+(z-z0)^2=validityRange
z=z(x)

the z0 in z(x) eliminate the one in the first equation so we dont need z0
in the arguments"""
    #by calcul
    a=(correction*alpha/(4*epsilon**2))**2
    b=-4*rf*a
    c=1+a*6*rf**2-correction**2*alpha*zf/(2*epsilon**2)
    d=-4*rf**3*a+correction**2*alpha*zf/(2*epsilon**2)*2*rf
    e=a*rf**4-correction**2*alpha*zf/(2*epsilon**2)*rf**2+correction**2*zf**2-validityRange**2
    return(a,b,c,d,e)

    
    


    

    
