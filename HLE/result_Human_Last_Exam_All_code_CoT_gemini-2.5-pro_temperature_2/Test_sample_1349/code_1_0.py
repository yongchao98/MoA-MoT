import numpy as np

# We found the supremum of X to be given by the expression 24 / (16 * pi^2 + 15).
# This script calculates the numerical value of this expression.

# The numbers in the final equation
numerator_val = 24
coefficient_of_pi_sq = 16
constant_term = 15

# Calculate pi squared
pi_squared = np.pi**2

# Calculate the denominator
denominator_val = coefficient_of_pi_sq * pi_squared + constant_term

# Calculate the supremum of X
sup_X = numerator_val / denominator_val

print("The supremum of X is given by the equation: numerator / (coefficient_of_pi_sq * pi^2 + constant_term)")
print(f"Numerator: {numerator_val}")
print(f"Coefficient of pi^2: {coefficient_of_pi_sq}")
print(f"Constant term: {constant_term}")
print(f"Using the value of pi^2 ~ {pi_squared:.4f}")
print(f"The denominator is {coefficient_of_pi_sq} * {pi_squared:.4f} + {constant_term} = {denominator_val:.4f}")
print(f"The numerical value for the supremum of X is approximately: {sup_X}")