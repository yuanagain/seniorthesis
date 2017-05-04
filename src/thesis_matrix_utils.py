## thesis_matrix_utils.py
## Yuan Wang

import numpy as np
from interval import *
from intervals.intervaln import *


def validateType(mat, show = False):
    """
    Validates data is properly formatted, esp to ensure that 
    intervals were not cast into lists or arrays
    """

    if type(mat) == np.matrix:
        if show == True:
            print(mat)
            print(type(mat))
            print(type(mat[0]))
            print(type(mat[0, 0]))

        if type(mat) != np.matrix:
            return False

        if type(mat[0]) != np.matrix:
            return False

        if type(mat[0, 0]) not in [interval, float, int]:
            return False

        return True

    elif type(mat) == np.array:
        if show == True:
            print(mat)
            print(type(mat))
            print(type(mat[0]))
            print(type(mat[0][0]))

        if type(mat) != np.matrix:
            return False

        if type(mat[0]) != np.matrix:
            return False

        if type(mat[0][0]) not in [interval, float, int]:
            return False

        return True

    if show == True:
        print(mat)
        print(type(mat))
        print(type(mat[0]))
        print(type(mat[0][0]))

    return False

def assertType(mat, show = True):
    assert(validateType(mat, show = show))

def midpoint(vect):
    """
    Takes midpoint
    """
    return toArray([[el[0].midpoint] for el in vect ])

def matrix_mult(mat, col):
    res = [0] * len(col[0])
    for i in range(len(mat)):
        row = mat[i]
        for j in range(len(j)):
            return

def contains(needle, haystack):
    intvln_1 = IntervalN(colToRow(needle))
    # print("Cols Cols Cols ")
    # print(self.__call__(Y_i) * h_k_int)
    intvln_2 = IntervalN(colToRow(haystack))
    # print("TYPES TYPES TYPES")
    # print(type(intvln_1))
    # print(type(intvln_2))
    return intvln_1 in intvln_2

def toArray(mat):
    """
    Converts list of lists or matrix to 2x2 numpy array
    """
    row_ct = len(mat)
    col_ct = len(mat[0])

    arr = np.empty(shape=(row_ct, col_ct), dtype=interval)

    for i in range(row_ct):
            for j in range(col_ct):
                if type(mat) == np.matrix:
                    arr[i, j] = mat[i, j]
                else:
                    arr[i, j] = mat[i][j]

    arr[:] = mat

    return arr

def toMatrix(mat):
    """
    Converts list of lists to 2x2 numpy array
    """
    row_ct = len(mat)
    col_ct = len(mat[0])

    arr = np.empty(shape=(row_ct, col_ct), dtype=interval)

    for i in range(row_ct):
            for j in range(col_ct):
                arr[i, j] = mat[i][j]

    arr[:] = mat

    return arr


def safeAdd(mat1, mat2):
    """
    because pyinterval and numpy don't play nice together. adds mat1 and mat2 pointwise
    """
    row_ct = len(mat1)
    col_ct = len(mat1[0])

    arr = np.empty(shape=(row_ct, col_ct), dtype=interval)

    for i in range(row_ct):
        for j in range(col_ct):
            # print(type(mat1[i][j]))
            # print(type(mat2[i][j]))
            # print(mat1[i][j])
            # print(mat2[i][j])
            arr[i][j] = mat1[i][j] + mat2[i][j]

    return arr


def matToArray(mat):
    """
    Converts list of lists to nxm numpy array
    """
    row_ct = len(mat)
    col_ct = len(mat[0])

    arr = np.empty(shape=(row_ct, col_ct), dtype=interval)

    for i in range(row_ct):
            for j in range(col_ct):
                arr[i, j] = mat[i][j]

    return arr

def listToCol(lst):
    """
    Converts list to column vector
    """

    return [[item] for item in lst]

def tupleToIntervalVector(lst):
    """
    Converts float tuple to column vector of intervals
    """
    return [[interval(item)] for item in lst]

def colToRow(lst):
    return [item[0] for item in lst]

    
            