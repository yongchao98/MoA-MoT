import numpy as np
from fractions import Fraction

# We are solving the system of linear equations Ax = b for the probabilities.
# The vector x represents the probabilities for the non-absorbing states:
# x = [p_ME, p_GH, p_CK, p_SP, p_L, p_KH, p_QR, p_TC]

# The coefficient matrix A is derived by rearranging the equations:
# e.g., p_ME = 0.5*p_GH + 0.5*p_CK  =>  1*p_ME - 0.5*p_GH - 0.5*p_CK = 0
A = np.array([
    [1, -1/2, -1/2, 0, 0, 0, 0, 0],               # p_ME equation
    [-1/4, 1, 0, -1/4, -1/4, -1/4, 0, 0],          # p_GH equation
    [-1/2, 0, 1, 0, 0, -1/2, 0, 0],               # p_CK equation
    [0, -1/4, 0, 1, 0, 0, 0, -1/2],               # p_SP equation (with p_TC substituted)
    [0, -1/3, 0, 0, 1, -1/3, -1/3, 0],             # p_L equation
    [0, -1/4, -1/4, 0, -1/4, 1, -1/4, 0],          # p_KH equation
    [0, 0, 0, 0, -1/3, -1/3, 1, 0],               # p_QR equation
    [0, 0, 0, -1/2, 0, 0, 0, 1]                    # p_TC equation
])

# The vector b contains the constant terms from the equations.
b = np.array([
    0,      # p_ME
    0,      # p_GH
    0,      # p_CK
    1/4,    # p_SP (from 1/4 * p_TR)
    0,      # p_L
    0,      # p_KH
    1/3,    # p_QR (from 1/3 * p_TR)
    0       # p_TC
])

# Solve the system of equations
probabilities = np.linalg.solve(A, b)

# Extract the probabilities for the relevant rooms
p_ME_float = probabilities[0]
p_GH_float = probabilities[1]
p_CK_float = probabilities[2]

# Convert floats to exact fractions
# The limit_denominator() function finds the simplest fraction that approximates the float
p_ME_frac = Fraction(p_ME_float).limit_denominator()
p_GH_frac = Fraction(p_GH_float).limit_denominator()
p_CK_frac = Fraction(p_CK_float).limit_denominator()

# Print the final calculation as requested
print("The probability of success from the Main Entrance (p_ME) is calculated based on the probabilities of its connected rooms:")
print("p_ME = (1/2) * p_GH + (1/2) * p_CK\n")
print("After solving the system, we find the probabilities:")
print(f"p_GH = {p_GH_frac}")
print(f"p_CK = {p_CK_frac}\n")
print("Substituting these values into the equation for p_ME:")
print(f"p_ME = (1/2) * ({p_GH_frac}) + (1/2) * ({p_CK_frac})")
print(f"p_ME = {1/2 * p_GH_frac} + {1/2 * p_CK_frac}")
print(f"p_ME = {p_ME_frac}\n")

print(f"The exact probability of reaching the Treasure Room from the Main Entrance is {p_ME_frac.numerator}/{p_ME_frac.denominator}.")

# Final answer in the required format
print(f'<<<{p_ME_frac.numerator}/{p_ME_frac.denominator}>>>')