## intervalmat.py
## Yuan Wang

import numpy as np

from intervaln import *
from __future__ import division
from interval import interval, inf, imath
import operator

class IntervalMat:
	def __init__(self, x = [[1,0,0,0],
							[0,1,0,0],
							[0,0,1,0],
							[0,0,0,1]] ):

		for i in range(len(x)):
			for j in range(len(x[0])):
				x[i][j] = interval(x[i][j])

		self.x = x

	def __mul__(self, other):
		if type(other) == IntervalN:

		else:
			row_ct = len(mat)
		    col_ct = len(mat[0])

		    arr = np.empty(shape=(row_ct, col_ct), dtype=interval)

		    for i in range(row_ct):
		            for j in range(col_ct):
		                arr[i, j] = mat[i][j]

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

    id = IntervalMat()

    print("============================")
    print("Construction")
    print("----------------------------")


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


if __name__=="__main__":
    main()