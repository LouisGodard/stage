from TraceARayRoots import *

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
    global distance
    distance=10000
    nbRay=21
    mini=int(-nbRay/2)
    maxi=-mini
    angleMax=np.pi/2
    tot=0
    for i in range (mini,maxi+1,1):
        rayon.direction=angleToVector(angleMax/maxi*i)
        abscisses,ordinates,compt=traceRay(distance)
        plt.plot(abscisses,ordinates,"k",linewidth=0.1)
        tot+=compt
        
        print(i,"and",angleMax/nbRay*i)
        
    print("the total number of peace of ray is :", tot)
    sol=[0,distance]
    plt.plot(sol,[0,0],"--")
    plt.ylim(0,1000)
    plt.xlim(0,distance)
    plt.show()
    plt.savefig("test1",transparent=True)
