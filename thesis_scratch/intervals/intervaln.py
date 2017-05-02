"""
Defines a basic n-dimensional interval and implements basic operations
on said intervals.

Author: 
Yuan Wang
Princeton University
"""

from __future__ import division
from interval import interval, inf, imath
import operator

class IntervalN(list): 
    def __init__(self, x = [interval([0, 0]), interval([0, 0]), interval([0, 0]), interval([0, 0])]):
        """
        Create defensive copy
        """
        self.x = [interval(intvl) for intvl in x]

    def clone(self):
        """
        Returns a clone of this interval
        """
        return IntervalN([interval(intvl) for intvl in self.x ])

    def midpoint(self):
        """
        Returns the midpoint of the hypercube
        """
        return IntervalN([intvl.midpoint for intvl in self.x ])

    def isZero(self):
        """
        Returns whether or not this interval is zero
        """
        zero = interval([0])
        is_zero = [intvl == zero for intvl in self.x]
        return sum(is_zero) == 0

    def __contains__(self, other):

        for i in range(len(self.x)):
            if other.x[i] not in self.x[i]:
                return False

        return True

    def hull(self):
        """
        Returns convex hull
        """
        return IntervalN( [intvl.hull([intvl]) for intvl in self.x] )

    def __add__(self, other):
        """
        Adds two hypercubes
        """
        return IntervalN( list(map(operator.sub, self.x, other.x)) )

    def __mul__(self, other):
        """
        Multiplies two hypercubes
        """
        if type(other) == IntervalN:
            raise TypeError("IntervalN cannot multiply with IntervalN") 

        return IntervalN( [intvl * other for intvl in self.x ] )


    def __truediv__(self, other):
        """
        Returns quotients of two hypercubes

        TODO: test more thoroughly
        """
        if type(other) == IntervalN:
            raise TypeError("IntervalN cannot multiply with IntervalN") 

        return IntervalN( [intvl / other for intvl in self.x ] )

        #return IntervalN( list(map(operator.truediv, self.x, other.x)) )

    def __sub__(self, other):
        """
        Subtracts other from self
        """
        return self.__add__(-other)

    def __neg__(self):
        """
        Returns negation of interval
        """
        return IntervalN([-intvl for intvl in self.x ])

    def __str__(self):
        """
        Returns string representation of complex interval
        """
        return "IntervalN(" + str(self.x) + ")"


    def __eq__(self, other):
        """
        Checks for equality between intervals
        """
        is_nequal = [self.x[i] != other.x[i] for i in range(len(self.x)) ]
        return sum(is_nequal) == 0

    def __ne__(self, other):
        """
        Checks for equality between intervals
        """
        return self.__eq__(other) == False

## Functions for generating copies of identities 

def _scalar(k = 0.0, dim = 4):
    """
    Returns a complex interval representing zero
    """
    return IntervalN( [interval([k]) for i in range(dim) ] )

def main():
    print("Testing IntervalN")
    xa = interval([1, 2])
    xb = interval([5, 6])
    xc = interval([-2, 4])
    xd = interval([-2, -9])

    x = IntervalN([xa, xb, xc, xd])
    x_clone = IntervalN([xa, xb, xc, xd])
    x_clone2 = x.clone()

    ya = interval([10, 2])
    yb = interval([52, 61])
    yc = interval([-22, 41])
    yd = interval([-34, -10])

    y = IntervalN([ya, yb, yc, yd])

    print("============================")
    print("Construction")
    print("----------------------------")

    print('x = ' + x.__str__())
    print('y = ' + y.__str__())

    print("============================")
    print("Comparisons")
    print("----------------------------")

    print(x_clone2 == x_clone)
    print(x_clone2 != x_clone)
    print(x_clone2 == y)
    print(x_clone2 != y)
    print(_scalar(7))

    print("============================")
    print("Operations")
    print("----------------------------")

    print(x + y)
    print(x - y)
    print(x * -2.1)
    print(x / 2)
    print(x.hull())


if __name__=="__main__":
    main()

