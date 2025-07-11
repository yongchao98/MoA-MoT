import numpy as np
from scipy.special import ai_zeros

# --- Step 1: Define constants and parameters ---
# Physical parameters from the problem
V0_eV = 15.0      # Potential parameter in eV
R_nm = 3.0        # Well radius in nm

# Fundamental constants in convenient units
hbar_c_eV_nm = 197.327  # Planck's constant * speed of light (eV*nm)
m_c2_eV = 511000.0      # Electron rest mass energy (eV)

# --- Step 2: Interpret and simplify the potential ---
# Based on the plan, we approximate the potential as a linear potential
# V(x) = c*x for x >= 0, where x = r - R.
# The slope 'c' is the derivative of V(r) at r=R.
# V(r) = V0 * (1 - (R/r)^2) for r>=R
# V'(r) = 2 * V0 * R^2 / r^3
# c = V'(R) = 2 * V0 / R
c_eV_per_nm = 2 * V0_eV / R_nm

# --- Step 3: Set up the energy level equation ---
# The energy levels for the "quantum bouncer" are given by:
# E_n = a_n * E_scale
# where -a_n are the zeros of the Airy function Ai(z), and E_scale is an energy factor.
# E_scale = (c^2 * hbar^2 / (2*m))^(1/3)
# We can calculate this using our convenient constants:
# E_scale = (c^2 * (hbar*c)^2 / (2*m*c^2))^(1/3)

# Calculate the energy scaling factor
term_inside_cube_root = (c_eV_per_nm**2 * hbar_c_eV_nm**2) / (2 * m_c2_eV)
E_scale_eV = term_inside_cube_root**(1.0/3.0)

# --- Step 4: Calculate the energy levels and their difference ---
# Find the first two zeros of the Airy function Ai(z) using SciPy.
# The function returns z_k where Ai(z_k) = 0. These are negative.
# Our a_n are the positive counterparts: a_n = -z_n.
num_zeros = 2
zeros_k = ai_zeros(num_zeros)[0]
a1 = -zeros_k[0]
a2 = -zeros_k[1]

# Calculate the first two energy levels
E1_eV = a1 * E_scale_eV
E2_eV = a2 * E_scale_eV

# Calculate the final energy difference
delta_E_eV = E2_eV - E1_eV

# --- Step 5: Print the final results ---
print("Based on a physical interpretation, the problem is solved as a 'quantum bouncer'.")
print(f"The potential is approximated by V(x) ≈ c * x, with c = 2*V₀/R.")
print(f"Slope c = 2 * {V0_eV} eV / {R_nm} nm = {c_eV_per_nm:.4f} eV/nm")
print(f"Energy scale factor = {E_scale_eV:.4f} eV")
print(f"Airy function zero a₁ = {a1:.4f}")
print(f"Airy function zero a₂ = {a2:.4f}")
print("\nCalculating the energy difference ΔE = E₂ - E₁:")
print(f"E₁ = {a1:.4f} * {E_scale_eV:.4f} eV = {E1_eV:.4f} eV")
print(f"E₂ = {a2:.4f} * {E_scale_eV:.4f} eV = {E2_eV:.4f} eV")
print("\nFinal Result:")
print(f"ΔE = {E2_eV:.4f} eV - {E1_eV:.4f} eV = {delta_E_eV:.4f} eV")

<<<2.7329>>>