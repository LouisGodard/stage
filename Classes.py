"""In this module, different classes are implemented.
The goal is to define different independant objects."""

class Source:
    """Source enables to create a ponctual source defined by its position
in the space (cartesian coordinates)"""

    def __init__(self, position=(0,0)):
        """position follow this format : (x,z)"""
        if type(position) is not tuple or len(position)!=2:
            raise FormatError("the position format is (x,z)")
        self.position=position
    

class Rayon:
    """A rayon object defines a ray by an origin and a direction at this
origin. The origin will be often at the same place of a source object"""

    def __init__(self, origin=(0,0), direction=(1,1)):
        if type(origin) is not tuple or len(origin)!=2:
            raise FormatError("The origin take the following form : (x,z)")
        if type(direction) is not tuple or len(direction)!=2 or\
           direction==(0,0): #we check that the direction is not the zero vector
            raise FormatError("""Direction has to be different
than the zero vector and is (x,z) coordinates""")
        self.origin=origin
        self.direction=direction


        
        
