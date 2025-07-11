import math

# Step 1: Define the constants from the problem's derivation.
# The solvability condition for the nonlinear boundary value problem leads to the
# following equation for the initial values (x₀, y₀, z₀) of the generating solution:
# 2 * (x₀² + y₀² + z₀²) * (1 - e⁻ᵀ) = α + β + γ
#
# Using the given constraint α + β + γ = 10²⁵ * (1 - e⁻ᵀ), we get:
# 2 * (x₀² + y₀² + z₀²) * (1 - e⁻ᵀ) = 10²⁵ * (1 - e⁻ᵀ)
#
# This simplifies to the equation of a sphere:
# x₀² + y₀² + z₀² = 10²⁵ / 2
# x₀² + y₀² + z₀² = 5 * 10²⁴
#
# This sphere is centered at the origin with a radius R, where R² = 5 * 10²⁴.

# Step 2: Calculate the area.
# The problem asks for the area of the surface defined by these values, which is
# the surface area of the sphere. The formula for the surface area of a sphere is
# A = 4 * π * R².

# The value of the radius squared
R_squared = 5 * 10**24

# The constants in the area formula
four = 4
pi = math.pi

# Calculate the surface area
surface_area = four * pi * R_squared

# Step 3: Print the final equation and the result.
print("The relationship between the initial values (x₀, y₀, z₀) defines a sphere.")
print(f"The equation of the sphere is: x₀² + y₀² + z₀² = {R_squared:.1e}")
print("\nThe area of this spherical surface is calculated as follows:")
print(f"Area = 4 * π * R²")
print(f"Area = {four} * {pi:.6f} * {R_squared:.1e}")
print(f"Area = {surface_area:.6e}")