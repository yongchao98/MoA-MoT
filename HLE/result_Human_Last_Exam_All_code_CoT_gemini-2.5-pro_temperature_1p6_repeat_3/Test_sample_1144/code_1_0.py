import math

# Define physical constants
h_bar = 1.054571817e-34  # Reduced Planck's constant (J·s)
a0 = 5.29177210903e-11   # Bohr radius (m)

# Given value
delta_x = 10e-12  # Uncertainty in position in meters (10 pm = 10 x 10^-12 m)

# 1. Calculate the uncertainty in momentum (delta_p)
# Using Heisenberg's Uncertainty Principle: delta_p = h_bar / (2 * delta_x)
delta_p = h_bar / (2 * delta_x)

# 2. Calculate the momentum of the electron in the first Bohr orbit (p)
# For the first Bohr orbit (n=1), p = h_bar / a0
p = h_bar / a0

# 3. Calculate the ratio of delta_p to p
ratio = delta_p / p

# Print the results
print(f"Given uncertainty in position (Δx): {delta_x:.2e} m")
print(f"Calculated uncertainty in momentum (Δp): {delta_p:.4e} kg·m/s")
print(f"Calculated momentum in the first Bohr orbit (p): {p:.4e} kg·m/s")
print("\nThe final ratio is calculated as follows:")
print(f"Ratio = Δp / p = {delta_p:.4e} / {p:.4e} = {ratio:.4f}")
print(f"\nThe ratio of the uncertainty of the momentum to the momentum is {ratio:.4f}.")