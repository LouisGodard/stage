from TraceARayPoly1d import *

def angleToVector(teta):
    """This function returns a vector which correspond to the same
direction than the angle (with x axis) given in argument"""
    y=0 #we are in 2D in the (x,z) plan
    x=1 #we fix x and we will determine z in consquence
    #we want z/x=tan(teta) so:
    z=np.tan(teta)*x
    return((x,y,z))
        

def rayShooting():
    """Shoots rays at different angles from the source.
A number of rays equals to nbRay will be shooted
between [-angleMax;angleMax]"""
    
    
    if nbRay==1:
        maxi=1
        mini=1
        peaceofAngle=angleMax
        #to trace one ray at angleMax
    else:
        maxi=(nbRay-1)/2
        mini=-maxi
        peaceofAngle=2*angleMax/(nbRay-1)
        #to trace rays at regular intervals between [-angleMax;angleMax] 

    tot=0 #to count the number of peace of ray
    indice=0 #to browse raysMatrix

    raysMatrix=np.empty(shape=(0,4),dtype=np.float64)
    
    for i in np.arange(mini,maxi+1,1):#maxi+1 to include maxi in the loop
        rayon=Rayon(source.position,angleToVector(peaceofAngle*i))
        ray,compt=traceRay(rayon)
        tot+=compt

        ray[0,0]=tot #like this the first indice of the ray which is not used
        #to describe a polynome, gives the last indice of the ray in raysMatrix
        
        raysMatrix=np.vstack((raysMatrix,ray))
        #the form of the ray matrix is a stack of peace of rays describe by
        #a,b,c,x1. the polynome of the peace of ray being ax^2+bx+c and the
        #abscisses of the limiting point being x1
        #when we meet a quadruple with a coefficient b or c infinite it means
        #we start a new ray
        
        
        indice+=1
        print(i,"and",peaceofAngle*i/np.pi,'π')
        
    print("the total number of peace of ray is :", tot)

    return(raysMatrix)

