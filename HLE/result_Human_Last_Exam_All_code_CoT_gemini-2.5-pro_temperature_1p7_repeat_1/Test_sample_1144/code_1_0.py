import scipy.constants as const

# Given uncertainty in position in meters
delta_x = 10e-12  # 10 pm

# Physical constants from scipy library
# hbar: reduced Planck constant
# a0: Bohr radius
hbar = const.hbar
bohr_radius = const.value('Bohr radius')

# Step 1: Calculate the uncertainty in momentum (delta_p)
# delta_p = hbar / (2 * delta_x)
delta_p = hbar / (2 * delta_x)

# Step 2: Calculate the momentum of the electron in the first Bohr orbit (p)
# p = hbar / a0 (where a0 is the radius of the first Bohr orbit)
p = hbar / bohr_radius

# Step 3: Calculate the ratio of delta_p to p
ratio = delta_p / p

# Also can be calculated directly as: ratio = bohr_radius / (2 * delta_x)
# This serves as a good check.

print(f"Uncertainty in position (Δx) = {delta_x:.2e} m")
print(f"Bohr radius (a₀) = {bohr_radius:.2e} m")
print(f"Reduced Planck constant (ħ) = {hbar:.2e} J·s")
print("-" * 30)
print(f"Uncertainty in momentum (Δp = ħ / (2·Δx)) = {delta_p:.4e} kg·m/s")
print(f"Momentum in first Bohr orbit (p = ħ / a₀) = {p:.4e} kg·m/s")
print("-" * 30)
print("The final equation is the ratio of these two values:")
print(f"{delta_p} / {p} = {ratio}")