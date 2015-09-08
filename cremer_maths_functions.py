def scalar_product(x, y):
	"""Perform scalar product on list x containing 6 values and list y containing 6 values to return value j"""
	z = []
    	for i in range(0,len(x)):
       		z.append( x[i] * y[i] )
	j = sum(z)
	return j
 
if __name__ == "__main__":
 
	a = [1, 2, 3, 4, 5, 6]

	b = [1, 2, 3, 4, 5, 6]

	c = scalar_product(a,b) 

	print(c)


def normalize(x):
	"""Performs a normalisation on the 6 values in the array ie sqrt(a^2 + b^2 + c^2 + d^2 + e^2 + f^2)""" 
	import math 
	squared = [y**2 for y in x] 
	sum_squared = sum(squared) 
	norm = math.sqrt(sum_squared)
	return norm 
	
