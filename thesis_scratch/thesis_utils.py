## thesis_utils.py
## Yuan Wang

import numdifftools as nd
import operator


def quad_distance(w, x, y, z):
    """mainly for plotting type 2"""
    return [w[i]**2 + x[i]**2 + y[i]**2 + z[i]**2 for i in range(len(w))]

def quad_sq_distance(x, y):
    """Computes the squared distance"""
    dists = [ x[i] - y[i] for i in range(len(x) )]
    dists = [ dists[i]**2 for i in range(len(x) )]
    return sum(dists)

def tuple_add(a, b):
    return tuple(map(operator.add, a, b) )
    
def tuple_subtract(a, b):
    b_neg = tuple([-k for k in b])
    return tuple(map(operator.add, a, b_neg) )
    
def list_subtract(a, b):
    return list(map(operator.sub, a, b))

def list_add(a, b):
    return list(map(operator.add, a, b))

def approx_derivs(x):
    """Approximate partial deritatives of x"""
    gx0 = g(*x)
    x_1_dot = ( g(*tuple_add(x, (dt, 0,  0 , 0 ) ) ) - gx0 ) / dt
    x_2_dot = ( g(*tuple_add(x, (0,  dt, 0 , 0 ) ) ) - gx0 ) / dt
    y_1_dot = ( g(*tuple_add(x, (0,  0,  dt, 0 ) ) ) - gx0 ) / dt
    y_2_dot = ( g(*tuple_add(x, (0,  0,  0 , dt) ) ) - gx0 ) / dt
    
    return (x_1_dot, x_2_dot, y_1_dot, y_2_dot)


def newton_iterate(x):
    gx0 = g(*x)
    x_1_dot = ( g(*tuple_add(x, (dt, 0,  0 , 0 ) ) ) - gx0 ) / dt
    x_2_dot = ( g(*tuple_add(x, (0,  dt, 0 , 0 ) ) ) - gx0 ) / dt
    y_1_dot = ( g(*tuple_add(x, (0,  0,  dt, 0 ) ) ) - gx0 ) / dt
    y_2_dot = ( g(*tuple_add(x, (0,  0,  0 , dt) ) ) - gx0 ) / dt

def g_nd(x):
    return g(*tuple(x))

def dots(x_0, lmbda):
    """
    dz1/dt = lambda_2 * z1^2 - (lambda_2 + lambda_3) * z1 * z2
    dz2/dt = lambda_1 * z2^2 - (lambda_1 + lambda_3) * z1 * z2
    http://www.math.kit.edu/iag3/~herrlich/seite/wws-11/media/wws-talk-valdez.pdf
    """
    x_1 = x_0[0]
    y_1 = x_0[1]
    x_2 = x_0[2]
    y_2 = x_0[3]

    # print(lmbda)
    lambda_1 = lmbda[0]
    lambda_2 = lmbda[1]
    lambda_3 = lmbda[2]

    x_1_dot = lambda_2 * (x_1**2 - y_1**2) - (lambda_2 + lambda_3) * (x_1*x_2 - y_1*y_2)
    y_1_dot = 2 * lambda_2 * x_1 * y_1 - (lambda_2 + lambda_3) * (x_1*y_2 + y_1*x_2)
    x_2_dot = lambda_1 * (x_2**2 - y_2**2) - (lambda_1 + lambda_3) * (x_1*x_2 - y_1*y_2)
    y_2_dot = 2 * lambda_1 * x_2 * y_2 - (lambda_1 +lambda_3) * (x_1*y_2 + y_1*x_2)

    return [x_1_dot, y_1_dot, x_2_dot, y_2_dot]