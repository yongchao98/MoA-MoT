import scipy.constants as const

# 1. Define constants and given values
# Reduced Planck constant (J·s)
hbar = const.hbar
# Bohr radius, which is the radius of the first Bohr orbit (m)
r1 = const.physical_constants['Bohr radius'][0]
# Given uncertainty in position (m)
delta_x = 10e-12  # 10 pm = 10 x 10^-12 m

# 2. Calculate the momentum of the electron in the first Bohr orbit (p)
# p = ħ / r₁
p = hbar / r1

# 3. Calculate the minimum uncertainty in momentum (Δp) using Heisenberg's Uncertainty Principle
# Δp * Δx ≥ ħ / 2  =>  Δp = ħ / (2 * Δx)
delta_p = hbar / (2 * delta_x)

# 4. Calculate the ratio of the uncertainty in momentum to the momentum
ratio = delta_p / p
# This can also be calculated with the simplified formula: ratio = r1 / (2 * delta_x)

# 5. Print the results
print("This script calculates the ratio of the uncertainty of an electron's momentum to its momentum in the first Bohr orbit.")
print("-" * 50)
print(f"Given uncertainty in position (Δx): {delta_x:.2e} m")
print(f"Radius of the first Bohr orbit (r₁): {r1:.4e} m")
print(f"Momentum of the electron (p = ħ/r₁): {p:.4e} kg·m/s")
print(f"Uncertainty in momentum (Δp = ħ/(2*Δx)): {delta_p:.4e} kg·m/s")
print("-" * 50)
print("The final ratio is calculated as Δp / p:")
print(f"Final Ratio = {delta_p:.4e} / {p:.4e} = {ratio:.4f}")

print(f"\nFinal Answer: {ratio}")
<<<2.6459>>>