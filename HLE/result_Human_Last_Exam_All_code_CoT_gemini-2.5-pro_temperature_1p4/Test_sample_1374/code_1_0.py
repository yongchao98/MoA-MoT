import math

# This script calculates the furthest distance from point A for a playdough shape
# that maximizes the gravitational field at A.

# The volume of the playdough in cubic meters is given.
V = 1.0

# The formula for the maximum distance (R_max) is derived from the shape's
# volume constraint. The optimal shape has a volume V = (4 * pi * R_max^3) / 15.
# Solving for R_max gives: R_max = (15 * V / (4 * pi))^(1/3).

# --- Calculation ---
numerator = 15 * V
denominator = 4 * math.pi
base = numerator / denominator
R_max = base**(1/3)

# --- Output ---
# As requested, we print the steps of the final calculation.
print("The equation for the furthest distance R_max is:")
print(f"R_max = (15 * V / (4 * π))^(1/3)")
print("\nSubstituting the given values:")
print(f"R_max = (15 * {V} / (4 * {math.pi}))^(1/3)")
print(f"R_max = ({numerator} / {denominator})^(1/3)")
print(f"R_max = ({base})^(1/3)")
print("\nThe furthest point on the surface of the playdough is:")
print(R_max)