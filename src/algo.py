import sys

max_iterates = 10
resolution = 10**-5

class rect:
	def __init__():
		this.perturb = 0 # track how much to perturb by if perturb bisect called.
		print('initialized')

	def bisect():
		"""
		Bisects the domain, returns two new rects.
		"""

	def perturb_bisect(perturb):
		"""
		Bisects the domain with the specified perturbation.
		"""

	def midpoint():
		"""
		Returns midpoint for rect
		"""


def integrate(rect, f):
	"""
	Integrate f over rect
	"""


def atResolution(rect):
	"""
	Returns true iff the rectangle described is at resolution
	"""

	return True

def recur(rect, f, zeroes):
	"""
	Populates zeroes with all zeroes of rect and f
	"""
	# Count zeroes
	numZeroes = argumentPrinciple(rect, f)

	# Nothing to add
	if (numZeroes == 0): return

	# Add to zeroes with multiplicity
	if (atResolution(rect)): 
		zeroes += (rect, numZeroes)
		return


	if (numZeroes == 1):
		if (newton(rect, f, zeroes)):





def newton(rect, f, zeroes):
	"""
	Searches for zero of f in rect, adds zero to zeroes.
	"""
	x = rect.midpoint()
	prev = rect.midpoint()
	converge = False
	iterates = 0
	while (converge == False):
		prev = x
		iterates += 1
		if iterates > max_iterates:
			return False

		# COMPUTE NEXT IN SEQUENCE HERE

		if prev == x: 
			converge = True

	numZeroes = argumentPrinciple(rect, f)


def argumentPrinciple(rect, f): 
	"""
	Computes N - P on the domain rect for f via argument principle
	"""
	return 0

def f(z):
	"""
	A complex function of the user's design.
	"""
	return 0

def main():
	zeroes = []
	recur(rect, f, zeroes)


if __name__=="__main__":
	main()