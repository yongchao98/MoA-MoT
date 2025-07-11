import math

# This script calculates the value of R based on the derived formula.

# The equation for the sphere of initial values is x₀² + y₀² + z₀² = R²,
# where R² = α / (1 - e⁻ᵀ).

# Given values and definitions:
# T = ln(10^34), so e^T = 10^34
# α = (1/2) * (e^(2T) - 1)

# To avoid numerical instability with very large and small numbers,
# we simplify the expression for R²:
# R² = [0.5 * (e^(2T) - 1)] / [1 - e⁻ᵀ]
# R² = [0.5 * (eᵀ - 1)(eᵀ + 1)] / [(eᵀ - 1)/eᵀ]
# R² = 0.5 * (eᵀ + 1) * eᵀ

# We use Python's arbitrary-precision integers for an exact calculation of R².
eT = 10**34

# Calculate R² using the stable, simplified formula.
# We use integer division // as the numerator is guaranteed to be even.
R_squared = (eT + 1) * eT // 2

# Calculate R by taking the square root.
R = math.sqrt(R_squared)

# As requested, we output the numbers in the final equation for the sphere.
# The equation is 1*x₀² + 1*y₀² + 1*z₀² = R²
print("The final equation for the set of initial values is:")
print(f"1 * x₀² + 1 * y₀² + 1 * z₀² = {R_squared}")
print("\nThe radius R of this sphere is:")
print(R)