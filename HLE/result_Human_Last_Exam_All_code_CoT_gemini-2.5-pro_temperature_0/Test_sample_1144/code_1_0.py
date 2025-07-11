import math

# Define the physical constants and given values
# Reduced Planck constant in J.s
hbar = 1.054571817e-34
# Bohr radius (radius of the first orbit, n=1) in meters
a0 = 5.29177210903e-11
# Uncertainty in position in meters (10 pm = 10 * 10^-12 m)
delta_x = 10e-12
# Principal quantum number for the first Bohr orbit
n = 1

# 1. Calculate the uncertainty in momentum (delta_p) using Heisenberg's Uncertainty Principle
# delta_p * delta_x >= hbar / 2
# We use the minimum uncertainty: delta_p = hbar / (2 * delta_x)
delta_p = hbar / (2 * delta_x)

# 2. Calculate the momentum (p) of the electron in the first Bohr orbit
# Angular momentum L = mvr = n * hbar. Since p = mv, then p = n * hbar / r.
# For the first orbit, n=1 and r = a0 (Bohr radius).
p = (n * hbar) / a0

# 3. Calculate the ratio of the uncertainty in momentum to the momentum
ratio = delta_p / p

# Print the results
print(f"Uncertainty in position (Δx): {delta_x:.2e} m")
print(f"Momentum in the first Bohr orbit (p): {p:.4e} kg·m/s")
print(f"Uncertainty in momentum (Δp): {delta_p:.4e} kg·m/s")
print("\nThe final equation is: Ratio = Δp / p")
print(f"Ratio = {delta_p:.4e} / {p:.4e}")
print(f"The calculated ratio is: {ratio}")

<<<2.6458860545150004>>>