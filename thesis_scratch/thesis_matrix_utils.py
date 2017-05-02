## thesis_matrix_utils.py
## Yuan Wang

import numpy as np
from interval import *

def matrix_mult(mat, col):
    res = [0] * len(col[0])
    for i in range(len(mat)):
        row = mat[i]
        for j in range(len(j)):
            return


def toArray(mat):
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
    row_ct = len(mat1)
    col_ct = len(mat1[0])

    arr = np.empty(shape=(row_ct, col_ct), dtype=interval)

    for i in range(row_ct):
        for j in range(col_ct):
            print(type(mat1[i][j]))
            print(type(mat2[i][j]))
            print(mat1[i][j])
            print(mat2[i][j])
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

    
            