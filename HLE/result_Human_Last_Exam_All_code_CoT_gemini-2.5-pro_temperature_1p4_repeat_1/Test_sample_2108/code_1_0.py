import math

# The problem is to find a specific ratio of radiated power for a rotating
# charged sphere. The physics of the problem leads to a purely geometric
# calculation based on the standard magnetic dipole radiation pattern.

# The final analytical expression for the maximum achievable ratio (R) is derived as:
# R = (pi / 6) * (16 - 7 * sqrt(2))

# This script calculates the numerical value of this expression.

# The numbers that constitute the final equation are:
# pi, 6, 16, 7, and 2.
pi_val = math.pi
term_divisor = 6
term_minuend = 16
term_subtrahend_coeff = 7
term_in_sqrt = 2

# Perform the calculation step-by-step
sqrt_val = math.sqrt(term_in_sqrt)
parenthesis_val = term_minuend - term_subtrahend_coeff * sqrt_val
final_result = (pi_val / term_divisor) * parenthesis_val

# Output the numbers used in the formula and the final result.
print("The final expression for the ratio is (pi / A) * (B - C * sqrt(D))")
print("The values for the components of this equation are:")
print(f"pi = {pi_val}")
print(f"A = {term_divisor}")
print(f"B = {term_minuend}")
print(f"C = {term_subtrahend_coeff}")
print(f"D = {term_in_sqrt}")
print("\nFinal calculated ratio:")
print(final_result)