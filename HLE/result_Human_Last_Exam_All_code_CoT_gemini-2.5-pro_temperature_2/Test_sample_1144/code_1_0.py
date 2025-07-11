import math

# Define the physical constants and given values.
# a0 is the Bohr radius in meters.
a0 = 5.29177e-11
# delta_x_pm is the given uncertainty in position in picometers.
delta_x_pm = 10

# Convert the uncertainty in position from picometers to meters for unit consistency.
delta_x_m = delta_x_pm * 1e-12

# As derived from the Heisenberg Uncertainty Principle and the Bohr model,
# the ratio of momentum uncertainty (Δp) to momentum (p) is:
# Ratio = a₀ / (2 * Δx)
ratio = a0 / (2 * delta_x_m)

# Print the explanation and the final equation with numerical values.
print("The ratio of the uncertainty in momentum (Δp) to the momentum (p) is calculated using the formula:")
print("Ratio = a₀ / (2 * Δx)\n")
print("Where:")
print(f"  a₀ (Bohr Radius) = {a0:.5e} m")
print(f"  Δx (Uncertainty in position) = {delta_x_pm} pm = {delta_x_m:.1e} m\n")

print("The final equation with these values is:")
print(f"Ratio = {a0:.5e} / (2 * {delta_x_m:.1e})")

# Print the final calculated result.
print("\nResult:")
print(f"The calculated ratio is: {ratio}")