import math

# Define the genera of the surfaces
g = 31
h = 17

# Step 1: Calculate the simplicial volume of the first surface, Sigma_g
# The formula for the simplicial volume of a surface with genus g >= 2 is ||Sigma_g|| = 4*g - 4.
sv_g = 4 * g - 4

# Step 2: Calculate the simplicial volume of the second surface, Sigma_h
# The formula for the simplicial volume of a surface with genus h >= 2 is ||Sigma_h|| = 4*h - 4.
sv_h = 4 * h - 4

# Step 3: Calculate the simplicial volume of the product Sigma_g x Sigma_h
# The formula is ||Sigma_g x Sigma_h|| = C(4, 2) * ||Sigma_g|| * ||Sigma_h||
# where C(4, 2) is the binomial coefficient "4 choose 2".
n = 2  # Dimension of the surfaces
m = 2
# Using math.comb for binomial coefficient, which is C(n+m, n)
coeff = math.comb(n + m, n)

# Calculate the final result
total_sv = coeff * sv_g * sv_h

# Print the step-by-step calculation
print(f"The simplicial volume of an oriented closed surface of genus g >= 2 is given by the formula ||\u03A3_g|| = 4g - 4.")
print(f"For g = {g}:")
print(f"||\u03A3_{g}|| = 4 * {g} - 4 = {sv_g}")
print("-" * 20)
print(f"For h = {h}:")
print(f"||\u03A3_{h}|| = 4 * {h} - 4 = {sv_h}")
print("-" * 20)
print(f"The simplicial volume of the product \u03A3_{g} x \u03A3_{h} is given by the formula:")
print(f"||\u03A3_{g} x \u03A3_{h}|| = C({n+m}, {n}) * ||\u03A3_{g}|| * ||\u03A3_{h}||")
print(f"where C({n+m}, {n}) = {coeff}.")
print("The final calculation is:")
print(f"||\u03A3_{g} x \u03A3_{h}|| = {coeff} * {sv_g} * {sv_h} = {total_sv}")