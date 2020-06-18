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
    global distance
    global delta
    global nbtot
    nbtot=400
    delta=0.1
    distance=200
    nbRay=2
    max=np.pi/13
    for i in range (nbRay):
        rayon.direction=angleToVector(max/nbRay*i)
        abscisses,ordinates=traceRay(distance)
        plt.plot(abscisses,ordinates,"k",linewidth=1)

        rayon.direction=angleToVector(-max/nbRay*i)
        abscisses,ordinates=traceRay(distance)
        plt.plot(abscisses,ordinates,"k",linewidth=1)

    sol=[0,distance]
    plt.plot(sol,[0,0],"--")
    plt.ylim(-2,10)
    plt.xlim(0,distance)
    plt.show()

