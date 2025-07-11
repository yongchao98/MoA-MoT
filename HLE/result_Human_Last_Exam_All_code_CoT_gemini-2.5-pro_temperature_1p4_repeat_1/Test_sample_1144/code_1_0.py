import math

# Step 1: Define the given values and physical constants.
# Uncertainty in position (Δx) in meters. 10 pm = 10 * 10^-12 m.
delta_x = 10e-12

# Bohr radius (a₀) in meters.
a0 = 5.29177e-11

# Step 2: The problem asks for the ratio of the uncertainty in momentum (Δp)
# to the momentum (p) of the electron in the first Bohr orbit.
# The ratio can be calculated using the simplified formula derived from
# the Uncertainty Principle and the Bohr model: Ratio = a₀ / (2 * Δx).

# Step 3: Calculate the ratio.
ratio = a0 / (2 * delta_x)

# Step 4: Print the final equation with the numbers substituted in.
print("To find the ratio of the uncertainty of the momentum to the momentum, we use the formula:")
print("Ratio = a₀ / (2 * Δx)\n")
print("Substituting the values:")
print(f"Ratio = {a0} m / (2 * {delta_x} m)")
print(f"Ratio = {a0} / {2 * delta_x}")
print(f"\nThe calculated ratio is: {ratio}")
