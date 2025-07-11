import numpy as np

# The problem requires finding the average of the complex coordinates z where the
# inverse of the gradient of S(z) approaches zero. This is equivalent to finding the
# coordinates z where the gradient of S(z) approaches infinity.
# These points correspond to the poles of the fields E(z) and B(z).
#
# The equation for B(z) is driven by a source term with a denominator:
# 4*z**4 - z**3 + z**2 + 1.
# The poles of B(z), and therefore the points z of interest, are the roots of the polynomial
# P(z) = 4*z**4 - z**3 + z**2 + 1 = 0.
#
# We need to find the average of these roots.
# For a polynomial of the form a_n*z**n + a_{n-1}*z**(n-1) + ... + a_0 = 0,
# the sum of the roots is given by Vieta's formulas as -a_{n-1}/a_n.
# The average of the roots is the sum divided by the degree n.

# Coefficients of the polynomial P(z) = 4*z**4 - 1*z**3 + 1*z**2 + 0*z + 1
a_4 = 4
a_3 = -1
a_2 = 1
a_1 = 0
a_0 = 1

# The degree of the polynomial
n = 4

# Calculate the sum of the roots using Vieta's formulas
sum_of_roots = -a_3 / a_4

# Calculate the average of the roots
average_of_roots = sum_of_roots / n

print("The coordinates z are the roots of the polynomial:")
print(f"{a_4}*z^4 + ({a_3})*z^3 + {a_2}*z^2 + {a_1}*z + {a_0} = 0\n")

print("The final calculation for the average value of z is:")
print(f"Average = (Sum of roots) / (Number of roots)")
print(f"Average = (-a_3 / a_4) / n")
print(f"Average = (-({a_3}) / {a_4}) / {n}")
print(f"Average = ({sum_of_roots}) / {n}")
print(f"Average = {average_of_roots}")

# For verification, we can also find the roots numerically using numpy and average them.
# This is not required for the final answer but serves as a good check.
# coeffs = [a_4, a_3, a_2, a_1, a_0]
# roots = np.roots(coeffs)
# numerical_average = np.mean(roots)
# print(f"\nNumerical verification: The average of the roots is {numerical_average}")
