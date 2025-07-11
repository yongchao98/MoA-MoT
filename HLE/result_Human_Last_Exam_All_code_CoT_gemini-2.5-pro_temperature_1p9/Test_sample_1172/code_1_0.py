import math

# --- User-defined Parameters ---
# These are example values. You can change them as needed.
# The condition d >> h must be satisfied for the approximation to be valid.

# Separation between wires in each circuit (meters)
h = 0.01
# Distance between the two circuits (meters)
d = 1.0
# Inner radius of the concentrator shells (meters)
R1 = 0.1
# Outer radius of the concentrator shells (meters)
R2 = 0.2

# --- Constants ---
# Permeability of free space (Henry per meter)
mu_0 = 4 * math.pi * 1e-7

# --- Calculations ---
# The formulas used are approximations valid for d >> h.

# M1': Mutual inductance per unit length for bare circuits
m1_prime = (mu_0 / (2 * math.pi)) * (h**2 / d**2)

# M2': Mutual inductance per unit length with concentrators
# The effective separation h is scaled by the factor (R2/R1).
m2_prime = (mu_0 / (2 * math.pi)) * (h**2 / d**2) * (R2 / R1)**2

# Delta_M': The change in mutual inductance per unit length (M1' - M2')
delta_m_prime = m1_prime - m2_prime

# --- Output ---
print("This script calculates the change in mutual inductance per unit length (Delta M')")
print("between two parallel wire pairs when surrounded by magnetic concentrators.")
print("The calculation is based on the approximation d >> h.\n")

print(f"Given parameters:")
print(f"  h  = {h} m")
print(f"  d  = {d} m")
print(f"  R1 = {R1} m")
print(f"  R2 = {R2} m\n")

print("The formula for the change in mutual inductance per unit length is:")
print("  ΔM' = (μ₀ * h² / (2 * π * d²)) * (1 - (R₂/R₁)²) \n")

print("--- Calculation Steps ---")

print(f"1. Mutual Inductance without concentrators (M₁'):")
print(f"   M₁' = (μ₀ / (2π)) * (h/d)²")
print(f"   M₁' = ({mu_0:.4e} / (2 * {math.pi:.4f})) * ({h}/{d})² = {m1_prime:.3e} H/m\n")

print(f"2. Mutual Inductance with concentrators (M₂'):")
print(f"   M₂' = M₁' * (R₂/R₁)²")
print(f"   M₂' = {m1_prime:.3e} * ({R2}/{R1})² = {m2_prime:.3e} H/m\n")

print(f"3. Change in Mutual Inductance (ΔM' = M₁' - M₂'):")
print(f"   Final equation with numbers:")
print(f"   ΔM' = ({mu_0:.4e} / (2 * {math.pi:.4f})) * ({h}² / {d}²) * (1 - ({R2}/{R1})²)")
print(f"   ΔM' = {m1_prime:.3e} H/m - {m2_prime:.3e} H/m")
print(f"   ΔM' = {delta_m_prime:.3e} H/m")

<<<(mu_0 * h**2 / (2 * pi * d**2)) * (1 - (R2/R1)**2)>>>