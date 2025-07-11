import math

# Let n be the number of objects in the category.
# The user can change this value to find the result for a different n.
n = 11

# The maximum number of eigenvalues greater than 2 is given by the formula floor((n+1)/3).
# We use integer division // which is equivalent to the floor function for positive results.
result = (n + 1) // 3

# Following the requested output format to print the final equation.
# e.g., for n=11, this will print "floor((11 + 1) / 3) = 4"
print(f"For n = {n}:")
print(f"floor(({n} + 1) / 3) = {result}")
