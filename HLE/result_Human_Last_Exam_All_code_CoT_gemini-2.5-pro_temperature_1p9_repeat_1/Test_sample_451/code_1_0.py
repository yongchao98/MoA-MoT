import math

# This script calculates the energy shifts due to centrifugal distortion for H-Br.
# The calculation follows these steps:
# 1. Define constants and input parameters in SI units.
# 2. Calculate the reduced mass (μ) of H-Br.
# 3. Calculate the moment of inertia (I).
# 4. Calculate the centrifugal distortion constant (D_J) in Joules.
# 5. Calculate the energy shift for each transition using ΔE = -4 * D_J * (J+1)³.
# 6. Convert the final energy shifts from Joules to quecto-electronvolts (qeV).

# --- Step 1 & 2: Define Constants and Inputs in SI units ---
r0 = 141.4e-12         # Bond length in meters
k = 400.0             # Force constant in N/m
m_H_amu = 1.008       # Mass of Hydrogen in amu
m_Br_amu = 79.904     # Mass of Bromine in amu

# Physical constants
AMU_TO_KG = 1.660539e-27
HBAR_Js = 1.054571817e-34
E_CHARGE_C = 1.602176634e-19  # Elementary charge in Coulombs for eV conversion

# --- Step 3: Calculate Reduced Mass (μ) ---
m_H_kg = m_H_amu * AMU_TO_KG
m_Br_kg = m_Br_amu * AMU_TO_KG
mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

# --- Step 4: Calculate Moment of Inertia (I) ---
I = mu_kg * r0**2

# --- Step 5: Calculate Centrifugal Distortion Constant (D_J) in Joules ---
# Formula: D_J = ħ⁴ / (2 * I² * k * r₀²)
D_J = (HBAR_Js**4) / (2 * I**2 * k * r0**2)

# --- Step 6 & 7: Calculate Shifts for Transitions and Convert to qeV ---
# Joules to qeV conversion factor: E(qeV) = E(J) * (1/e) * 10³⁰
JOULES_TO_QEV_FACTOR = 1e30 / E_CHARGE_C

# --- Calculation for the J=0 to J=1 transition ---
J1 = 0
shift_j1_joules = -4 * D_J * (J1 + 1)**3
shift_j1_qev = shift_j1_joules * JOULES_TO_QEV_FACTOR

# --- Calculation for the J=1 to J=2 transition ---
J2 = 1
shift_j2_joules = -4 * D_J * (J2 + 1)**3
shift_j2_qev = shift_j2_joules * JOULES_TO_QEV_FACTOR

# --- Step 8: Print the detailed results ---
print("--- Intermediate Values ---")
print(f"Reduced mass (μ): {mu_kg:.4e} kg")
print(f"Moment of inertia (I): {I:.4e} kg·m²")
print(f"Centrifugal distortion constant (D_J): {D_J:.4e} J")
print(f"Conversion factor from Joules to qeV: {JOULES_TO_QEV_FACTOR:.4e} qeV/J")

print("\n--- Energy Shift for J=0 to J=1 transition ---")
print(f"The calculation is: ΔE_shift = -4 * D_J * (J+1)³")
print(f"ΔE_shift (in Joules) = -4 * {D_J:.4e} J * ({J1}+1)³ = {shift_j1_joules:.4e} J")
print(f"ΔE_shift (in qeV) = {shift_j1_joules:.4e} J * {JOULES_TO_QEV_FACTOR:.4e} qeV/J")
print(f"Final Answer 1: The energy shift is {shift_j1_qev:.4e} qeV.")


print("\n--- Energy Shift for J=1 to J=2 transition ---")
print(f"The calculation is: ΔE_shift = -4 * D_J * (J+1)³")
print(f"ΔE_shift (in Joules) = -4 * {D_J:.4e} J * ({J2}+1)³ = {shift_j2_joules:.4e} J")
print(f"ΔE_shift (in qeV) = {shift_j2_joules:.4e} J * {JOULES_TO_QEV_FACTOR:.4e} qeV/J")
print(f"Final Answer 2: The energy shift is {shift_j2_qev:.4e} qeV.")
