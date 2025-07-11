import numpy as np

# The value of the integral is given by the expression I = pi * ln(1 + sqrt(2)).
# This script calculates the numerical value of this expression and shows the intermediate components.

# Define the constants from the formula
pi = np.pi
sqrt2 = np.sqrt(2)

# Calculate the term inside the natural logarithm
term_in_ln = 1 + sqrt2

# Calculate the natural logarithm
ln_value = np.log(term_in_ln)

# Calculate the final value of the integral I
integral_value = pi * ln_value

# As requested, here are the values of each number in the final equation.
print("The analytical solution for the integral is I = pi * ln(1 + sqrt(2)).")
print("Let's break down the final calculation:")
print(f"pi ≈ {pi}")
print(f"sqrt(2) ≈ {sqrt2}")
print(f"1 + sqrt(2) ≈ {term_in_ln}")
print(f"ln(1 + sqrt(2)) ≈ {ln_value}")
print("\nFinal calculation:")
print(f"I = {pi} * {ln_value}")
print(f"I ≈ {integral_value}")