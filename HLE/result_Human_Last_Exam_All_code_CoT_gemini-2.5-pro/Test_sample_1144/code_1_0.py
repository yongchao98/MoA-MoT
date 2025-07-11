import math

# Define the physical constants
# Reduced Planck constant in J·s
h_bar = 1.054571817e-34
# Bohr radius in meters
a0 = 5.291772109e-11

# Given uncertainty in position in meters (10 pm = 10 * 10^-12 m)
delta_x = 10e-12

# Step 1: Calculate the uncertainty in momentum (Δp)
delta_p = h_bar / (2 * delta_x)

# Step 2: Calculate the momentum of the electron in the first Bohr orbit (p)
p = h_bar / a0

# Step 3: Calculate the ratio
ratio = delta_p / p

# Print the results
print("This script calculates the ratio of the momentum uncertainty to the momentum for an electron in the first Bohr orbit.")
print(f"Given position uncertainty (Δx): {delta_x} m")
print(f"Bohr radius (a₀): {a0} m")
print("-" * 30)
print(f"Calculated momentum uncertainty (Δp): {delta_p:.4g} kg·m/s")
print(f"Calculated momentum in the first Bohr orbit (p): {p:.4g} kg·m/s")
print("-" * 30)
print("The final ratio is Δp / p:")
print(f"{delta_p:.4g} / {p:.4g} = {ratio}")