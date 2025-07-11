from fractions import Fraction

# Parameters from the problem
d = 2  # Dimension of the ambient space R^d
s = Fraction(8, 5)  # Dimension of the Frostman measure

# Perform the calculation step-by-step based on Wolff's theorem.
# The exponent for the decay of the L^2 norm squared (the integral) is -(d-s).
# The exponent for the L^2 norm itself is half of that. This is c.

# Step 1: Calculate d - s
step1_val = d - s
# Step 2: Calculate the exponent for the integral, -(d-s)
integral_exponent = -step1_val
# Step 3: Calculate the exponent for the L^2 norm, c
c = integral_exponent / 2

# Print the final equation, showing all the numbers used in the calculation.
print(f"c = -(d - s) / 2 = -({d} - {s}) / 2 = -({step1_val}) / 2 = {integral_exponent} / 2 = {c}")