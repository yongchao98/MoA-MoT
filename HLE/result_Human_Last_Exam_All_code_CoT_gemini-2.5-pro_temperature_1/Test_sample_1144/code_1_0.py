import math

# Step 1: Define the given values and physical constants.
# Uncertainty in position in meters (10 pm = 10e-12 m)
delta_x = 10e-12

# Bohr radius in meters
a0 = 5.29e-11

print("This script calculates the ratio of the uncertainty of momentum to the momentum of an electron in the first Bohr orbit.")
print("-" * 80)

# Step 2: Use the simplified formula derived from Heisenberg's Uncertainty Principle and the Bohr model.
# The ratio (Δp / p) simplifies to a₀ / (2 * Δx).
# Δp = ħ / (2 * Δx)
# p = ħ / a₀
# Δp / p = (ħ / (2 * Δx)) / (ħ / a₀) = a₀ / (2 * Δx)

print(f"Given uncertainty in position (Δx) = {delta_x:.2e} m")
print(f"The Bohr radius (a₀) for the first orbit is {a0:.2e} m")
print("\nThe ratio is calculated using the formula: Ratio = a₀ / (2 * Δx)")
print("-" * 80)

# Step 3: Calculate the ratio.
ratio = a0 / (2 * delta_x)

# Step 4: Print the final equation with the values and the result.
print("Plugging in the values into the equation:")
print(f"Ratio = {a0} / (2 * {delta_x})")
print(f"Ratio = {ratio}")

<<<2.645>>>