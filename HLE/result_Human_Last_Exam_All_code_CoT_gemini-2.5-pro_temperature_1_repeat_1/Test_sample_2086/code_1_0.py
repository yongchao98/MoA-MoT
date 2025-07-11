import math

# Let n be the number of objects. We set a specific value for n to demonstrate the calculation.
# This value can be changed by the user.
n = 8

# Based on the analysis, the maximum number of eigenvalues greater than 2
# for the matrix M(C_n) is given by the formula floor((n + 1) / 3).
# This formula is derived by finding the optimal partition of n objects
# into isomorphism classes (using only classes of size 2 and 3)
# that maximizes the count of eigenvalues > 2.

# We apply the formula. In Python, integer division // performs the floor operation.
numerator = n + 1
denominator = 3
max_eigenvalues = numerator // denominator

# The final output is printed, showing the formula with the specific numbers used.
print(f"For a given n, the maximum number of eigenvalues greater than 2 is calculated using the formula:")
print(f"result = floor((n + 1) / 3)")
print(f"For n = {n}, the calculation is:")
print(f"result = floor(({n} + 1) / {denominator}) = floor({numerator} / {denominator}) = {max_eigenvalues}")