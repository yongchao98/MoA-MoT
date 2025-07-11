import scipy.constants as const

# --- Given values and constants ---

# Uncertainty in position (Δx) is 10 pm. Convert pm to meters.
# 1 pm = 1e-12 m
delta_x_pm = 10
delta_x = delta_x_pm * 1e-12  # in meters

# Reduced Planck constant (ħ) in J·s
h_bar = const.hbar

# Radius of the first Bohr orbit (r₁) is the Bohr radius (a₀) in meters
r1 = const.physical_constants['Bohr radius'][0]

# --- Step 1: Calculate the uncertainty in momentum (Δp) ---
# From Heisenberg Uncertainty Principle: Δp = ħ / (2 * Δx)
delta_p = h_bar / (2 * delta_x)

# --- Step 2: Calculate the momentum in the first Bohr orbit (p) ---
# From Bohr model: p = ħ / r₁
p = h_bar / r1

# --- Step 3: Calculate the ratio of the uncertainties ---
ratio = delta_p / p

# --- Step 4: Print the results step-by-step ---
print(f"Step 1: Calculate the uncertainty in momentum (Δp)")
print(f"Given uncertainty in position (Δx) = {delta_x_pm} pm = {delta_x:.2e} m")
print(f"Using the Heisenberg Uncertainty Principle: Δp = ħ / (2 * Δx)")
print(f"Δp = {h_bar:.4e} J·s / (2 * {delta_x:.2e} m)")
print(f"Δp = {delta_p:.4e} kg·m/s\n")

print(f"Step 2: Calculate the momentum in the first Bohr orbit (p)")
print(f"Radius of the first Bohr orbit (r₁) = {r1:.4e} m")
print(f"Using the Bohr model: p = ħ / r₁")
print(f"p = {h_bar:.4e} J·s / {r1:.4e} m")
print(f"p = {p:.4e} kg·m/s\n")

print(f"Step 3: Calculate the final ratio (Δp / p)")
print(f"Ratio = Δp / p = {delta_p:.4e} / {p:.4e}")
# The final result can also be simplified to r₁ / (2 * Δx)
print(f"Ratio = r₁ / (2 * Δx) = {r1:.4e} m / (2 * {delta_x:.2e} m)")
print(f"The final calculated ratio is: {ratio}")

<<<2.645886036838896>>>