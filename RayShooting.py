from TraceARayPoly1d import *


def angleToVector(teta):
    """This function returns a vector which correspond to the same
direction than the angle (with x axis) given in argument"""
    x=1 #we fix x and we will determine z in consquence
    #we want z/x=tan(teta) so:
    z=np.tan(teta)*x
    return((x,z))
        

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
    indice=0 #to browse raysIndex

    raysMatrix=np.empty(shape=(0,5),dtype=np.float64)#will contain all the rays in a row
    raysIndex=np.empty(shape=(nbRay,),dtype=np.int16)#indexation of the rays in raysMatrix
    
    for i in np.arange(mini,maxi+1,1):#put maxi+1 to include maxi in the loop
        
        rayon=Rayon(source.position,angleToVector(peaceofAngle*i))#rayon is
        #the ray we will trace
        ray,compt=traceRay(rayon)
        tot+=(compt+1)

        
        raysIndex[indice]=tot #the rays index contains the indice just above
        #of the end of the i th ray

        raysMatrix=np.vstack((raysMatrix,ray))
        #the form of the ray matrix is a stack of peace of rays describe by
        #a,b,c,x1,reflexion. the polynome of the peace of ray being ax^2+bx+c and the
        #abscisses of the limiting point being x1, reflexion indicating if a reflexion happened
        #when we meet a 5-uple with a coefficient b or c infinite it means
        #a new ray begin
        
        indice+=1
        print("ray at indice",i,"and at angle",peaceofAngle*i/np.pi*180,'degree(s)')
        
    print("the total number of peaces of ray is :", tot)

    return(raysMatrix,raysIndex)

