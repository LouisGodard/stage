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
    """Shoots rays at different angles at the source.
A number of rays equals to nbRay will be shooted
between [-angleMax;angleMax]
if nbRay is not an odd number you will have one ray more than expected
(the one with the horizontal angle)"""
    
    
    if nbRay==1:
        maxi=1
        mini=1
        peaceofAngle=angleMax
        #to trace one ray at angleMax
    else:
        maxi=(nbRay-1)/2
        mini=-maxi
        peaceofAngle=2*angleMax/(nbRay-1)

    tot=0 #to count the number of peace of ray
    indice=0 #to browse in raysMAtrix
    raysMatrix=np.empty(shape=(2,int(maxi+1-mini)),dtype=object)
    
    for i in np.arange(mini,maxi+1,1):#maxi+1 to include maxi in the loop
        rayon.direction=angleToVector(peaceofAngle*i)
        abscisses,ordinates,compt=traceRay(distance)
        tot+=compt
        raysMatrix[0,indice]=abscisses
        raysMatrix[1,indice]=ordinates
        indice+=1
        print(i,"and",peaceofAngle*i/np.pi,'π')
        
    print("the total number of peace of ray is :", tot)

    return(raysMatrix)

