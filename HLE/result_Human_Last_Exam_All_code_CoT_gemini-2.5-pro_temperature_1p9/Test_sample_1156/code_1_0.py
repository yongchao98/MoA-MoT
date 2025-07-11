import numpy as np

# This script calculates and prints the normalised density of the invariant measure.
# The formula is derived using the conjugacy of the map T(x) to the Gauss map.
# The final density is given by rho(x) = 1 / (C * (x^p1 + x^p2)),
# where C is a constant derived from ln(2).

# 1. Define the parameters of the final equation
numerator = 1.0
power1 = 0.5
power2 = 1.0
constant_factor = 2.0

# 2. Calculate the constant C = 2 * ln(2)
C = constant_factor * np.log(2)

# 3. Print the results in a clear format.
# The code outputs each number in the final equation as requested.
print("The normalised density of the invariant measure is rho(x).")
print("The formula for the density is:")
print(f"rho(x) = {numerator} / (C * (x^{power1} + x^{power2}))")
print("\nWhere the constant C is calculated as:")
print(f"C = {constant_factor} * ln(2) = {C}")
print("\nSo the final equation for the density is:")
print(f"rho(x) = {numerator} / ({C} * (x**{power1} + x**{power2}))")