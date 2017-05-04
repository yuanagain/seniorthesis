## systems.py

def evolution_1(x_0, lmbda):
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


class evolution_1a:
	def __init__(self, lmbda)
		self.lmbda = lmbda

	def __call__(self, x_0):

		x_1_dot = lambda_2 * (x_1**2 - y_1**2) - (lambda_2 + lambda_3) * (x_1*x_2 - y_1*y_2)
    	y_1_dot = 2 * lambda_2 * x_1 * y_1 - (lambda_2 + lambda_3) * (x_1*y_2 + y_1*x_2)
	    x_2_dot = lambda_1 * (x_2**2 - y_2**2) - (lambda_1 + lambda_3) * (x_1*x_2 - y_1*y_2)
    	y_2_dot = 2 * lambda_1 * x_2 * y_2 - (lambda_1 +lambda_3) * (x_1*y_2 + y_1*x_2)

		return return [x_1_dot, y_1_dot, x_2_dot, y_2_dot]

	def __str__(self):
		return 	"dz1/dt = lambda_2 * z1^2 - (lambda_2 + lambda_3) * z1 * z2 \n" + \
				"dz2/dt = lambda_1 * z2^2 - (lambda_1 + lambda_3) * z1 * z2" + \ 
				"lambda_1: " + str(self.lambda_1) + \
				"; lambda_2: " + str(self.lambda_2) + \
				"; lambda_3: " + str(self.lambda_3)