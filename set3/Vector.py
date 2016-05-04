# Vector.py
# Aritra Biswas

import math

class Vector(list):
    '''A list-based vector class capable of element-wise addition and
    subtraction of two vectors, element-wise addition and subtraction of
    constants, i.e. c + v is interpreted as (c, c, ...) + v, and
    multiplication or division by a constant.'''

    def __add__(self, other):
        '''Add either another vector (element-by-element) or a constant,
        equivalent to adding the vector constant*[1, 1, ..., 1].'''

        try:
            return Vector(map(lambda x, y: x + y, self, other))
        except TypeError:
            return Vector(map(lambda x: x + other, self))

    def __radd__(self, other):
        '''Implement c + vec as vec + c.'''

        return self + other
    
    def __sub__(self, other):
        '''Subtract either another vector (element-by-element) or a constant,
        equivalent to subtracting the vector constant*[1, 1, ..., 1].'''

        try:
            return Vector(map(lambda x, y: x - y, self, other))
        except TypeError:
            return Vector(map(lambda x: x - other, self))

    def __rsub__(self, other):
        '''Implement c - vec as -(vec - c).'''

        return -1 * (self - other)

    def __mul__(self, c):
        '''Multiply by a constant.'''

        return Vector([c * x for x in self])

    def __rmul__(self, c):
        '''Implement c * vec as vec * c.'''

        return self * c
    
    def __div__(self, c):
        '''Divide by a constant.'''

        return Vector([x / c for x in self])

    def __abs__(self):
        '''Get the magnitude.'''

        return math.sqrt(sum([x*x for x in self]))
