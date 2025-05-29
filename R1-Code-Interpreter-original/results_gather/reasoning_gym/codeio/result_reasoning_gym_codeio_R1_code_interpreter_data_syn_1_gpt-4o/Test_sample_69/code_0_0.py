import math
from fractions import Fraction

# Given output dimensions
dimensions = [6, 8, 1, 9, 4]

# Calculate the product of the dimensions
product = math.prod(dimensions)

# Derive the ratios by dividing each dimension by the greatest common divisor
gcd = math.gcd(*dimensions)
ratios = [Fraction(d, gcd) for d in dimensions]

# Convert ratios to a string
rstring = ":".join(str(r.numerator) for r in ratios)

# Calculate a feasible volume
volume = product / gcd**len(dimensions)

print(product, rstring, volume)