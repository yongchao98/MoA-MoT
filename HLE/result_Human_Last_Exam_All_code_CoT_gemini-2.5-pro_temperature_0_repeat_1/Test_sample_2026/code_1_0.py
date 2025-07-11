import math

# From the derivation, the initial values (x₀, y₀, z₀) must satisfy the equation
# of a sphere: x₀² + y₀² + z₀² = R², where R² = 5 * 10²⁴.
# The problem asks for the area of this surface.

# The radius squared of the sphere.
R_squared = 5 * 10**24

# The formula for the surface area of a sphere is A = 4 * π * R².
# We can write the final equation for the area as A = 20 * π * 10²⁴.
coeff = 20
power = 24

# Calculate the final area.
surface_area = coeff * math.pi * (10**power)

print("The condition on the initial values (x₀, y₀, z₀) defines a sphere with the equation:")
print(f"x₀² + y₀² + z₀² = {R_squared:.0e}")
print("\nThe area of this sphere is calculated using the formula A = 4 * π * R².")
print(f"A = 4 * π * ({R_squared:.0e})")
print(f"A = {coeff} * π * 10^{power}")
print(f"\nThe final calculated area is: {surface_area}")