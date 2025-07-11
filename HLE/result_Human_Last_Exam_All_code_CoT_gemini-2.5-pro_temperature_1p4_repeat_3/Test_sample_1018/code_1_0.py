import math

# Step 1: Define the genera of the two surfaces.
g = 31
h = 17

# Step 2: Compute the simplicial volume for each surface.
# The formula for the simplicial volume of a surface of genus g >= 1 is ||Σ_g|| = 4g - 4.
sv_g = 4 * g - 4
sv_h = 4 * h - 4

print(f"The simplicial volume of the first surface, ||Σ_{g}||, is:")
print(f"4 * {g} - 4 = {sv_g}")

print(f"\nThe simplicial volume of the second surface, ||Σ_{h}||, is:")
print(f"4 * {h} - 4 = {sv_h}")

# Step 3: Compute the simplicial volume of the product.
# The formula is ||Σ_g x Σ_h|| = C(4, 2) * ||Σ_g|| * ||Σ_h||.
# The binomial coefficient C(4, 2) is 6.
coeff = math.comb(4, 2)
total_sv = coeff * sv_g * sv_h

print(f"\nThe simplicial volume of the product, ||Σ_{g} x Σ_{h}||, is:")
print(f"{coeff} * {sv_g} * {sv_h} = {total_sv}")