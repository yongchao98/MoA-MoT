import numpy as np

# The problem is to find the maximum real number c that satisfies the given inequality
# for all positive integers n. This value corresponds to the infimum of the
# minimum eigenvalues of the matrices A_n, where (A_n)_ij = n - |i-j|.
#
# The infimum is the limit of the minimum eigenvalues as n approaches infinity.
# This limit is known to be 4 / pi^2.
#
# This script calculates the numerical value of c from this expression.
# The final equation is c = 4 / pi^2.

# Define the numbers in the equation
numerator = 4
base = np.pi
exponent = 2

# Calculate the value of c
c = numerator / (base ** exponent)

# Per the instruction, printing the numbers involved in the final equation
print(f"The expression for c is: {numerator} / ({base}^'pi'^{exponent})")
print(f"The numerator is: {numerator}")
print(f"The base in the denominator is pi, which is approximately: {base}")
print(f"The exponent in the denominator is: {exponent}")
print(f"The calculated value of c is: {c}")
