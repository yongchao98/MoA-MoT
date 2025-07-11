import math

# Define the constants and given values in picometers (pm)
# Bohr radius (a0)
a0_pm = 52.9177
# Uncertainty in position (delta_x)
delta_x_pm = 10.0

# The ratio is calculated using the simplified formula: a0 / (2 * delta_x)
ratio = a0_pm / (2 * delta_x_pm)

# Print the explanation and the final calculation
print("To find the ratio of the uncertainty of the momentum (Δp) to the momentum (p) in the first Bohr orbit, we use the formula:")
print("Ratio = Δp / p = a₀ / (2 * Δx)")
print("\nWhere:")
print(f"  a₀ (Bohr radius) = {a0_pm} pm")
print(f"  Δx (Uncertainty in position) = {delta_x_pm} pm")
print("\nPlugging in the values, the final equation is:")
print(f"Ratio = {a0_pm} / (2 * {delta_x_pm})")
print(f"\nResult:")
print(f"The calculated ratio is: {ratio}")