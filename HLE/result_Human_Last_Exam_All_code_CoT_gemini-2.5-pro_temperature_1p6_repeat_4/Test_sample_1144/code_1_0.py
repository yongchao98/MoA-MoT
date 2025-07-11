import math

# Step 1: Define physical constants and given values
hbar = 1.054e-34  # Reduced Planck's constant in J*s
bohr_radius = 5.29e-11  # Bohr radius (radius of the first orbit) in meters
delta_x = 10e-12  # Uncertainty in position in meters (10 pm)

# Step 2: Calculate the uncertainty in momentum (delta_p)
# Using Heisenberg's Uncertainty Principle: delta_p = hbar / (2 * delta_x)
delta_p = hbar / (2 * delta_x)

# Step 3: Calculate the momentum of the electron in the first Bohr orbit (p1)
# Using the Bohr model: p1 = hbar / bohr_radius
p1 = hbar / bohr_radius

# Step 4: Calculate the ratio of the uncertainty in momentum to the momentum
ratio = delta_p / p1

# Print the values used in the final ratio calculation and the result
print(f"The uncertainty in momentum (Δp) is: {delta_p:.4g} kg*m/s")
print(f"The momentum in the first Bohr orbit (p1) is: {p1:.4g} kg*m/s")
print(f"The final ratio (Δp / p1) is calculated as: {delta_p:.4g} / {p1:.4g}")
print(f"Resulting Ratio = {ratio}")
