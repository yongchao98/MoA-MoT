import numpy as np

# The equation for the generating amplitude c1 is derived from the bifurcation analysis.
# The equation is: 2*pi - (c1**2 / 2) * (exp(4*pi) - 1) = 0
# We can write this in the form: A - B * c1**2 = 0

# Calculate the coefficients A and B
A = 2 * np.pi
B = (np.exp(4 * np.pi) - 1) / 2

# The problem requires outputting each number in the final equation.
print("The equation for the generating amplitude c1 is of the form A - B * c1**2 = 0:")
print(f"A = {A}")
print(f"B = {B}")
print(f"So the equation is: {A} - {B} * c1**2 = 0")

# Solve for c1. We are looking for the positive root c1 > 0.
# c1**2 = A / B
c1_squared = A / B
c1 = np.sqrt(c1_squared)

print("\nThe value of the first positive root c1 is:")
print(c1)