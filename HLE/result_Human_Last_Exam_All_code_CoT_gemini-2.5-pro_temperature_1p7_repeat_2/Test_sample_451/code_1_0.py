import math

# --- Physical and Conversion Constants ---
# Given values
r0_m = 141.4e-12         # Bond length in meters
k_Npm = 400.0             # Force constant in N/m
m_H_amu = 1.008           # Atomic mass of Hydrogen in amu
m_Br_amu = 79.904         # Atomic mass of Bromine in amu

# Fundamental constants and conversion factors
amu_to_kg = 1.660539e-27  # kg/amu
hbar_Js = 1.0545718e-34   # Reduced Planck constant in J·s
e_charge_C = 1.6021766e-19 # Elementary charge in Coulombs
J_to_eV = 1 / e_charge_C  # Conversion factor from Joules to eV
eV_to_qeV = 1e30          # Conversion factor from eV to qeV

# --- Step 1: Calculate Reduced Mass (μ) ---
m_H_kg = m_H_amu * amu_to_kg
m_Br_kg = m_Br_amu * amu_to_kg
mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

# --- Step 2: Calculate the Centrifugal Distortion Constant (D_E) in Joules ---
# First, we need to find the moment of inertia (I), the rotational constant (B_E),
# and the vibrational angular frequency (ω).

# Moment of inertia, I = μ * r₀²
I = mu_kg * r0_m**2

# Rotational constant in Joules, B_E = ħ² / (2 * I)
B_E = hbar_Js**2 / (2 * I)

# Vibrational angular frequency in rad/s, ω = sqrt(k / μ)
omega = math.sqrt(k_Npm / mu_kg)

# Centrifugal distortion constant in Joules, D_E = 4 * B_E³ / (ħω)²
D_E = (4 * B_E**3) / (hbar_Js**2 * omega**2)

# --- Step 3: Calculate and Print Energy Shifts ---
# The formula for the energy shift is ΔE_shift = -4 * D_E * (J+1)³

print("The formula for the energy shift of a rotational transition J -> J+1 is:")
print("ΔE_shift = -4 * D_E * (J+1)³\n")
print(f"Calculated values:")
print(f"Reduced mass, μ = {mu_kg:.4e} kg")
print(f"Centrifugal distortion constant, D_E = {D_E:.4e} J\n")

# --- Transition 1: J = 0 to J = 1 ---
J1 = 0
j_term_1 = (J1 + 1)**3
shift_1_J = -4 * D_E * j_term_1
shift_1_qeV = shift_1_J * J_to_eV * eV_to_qeV

print("1. For the transition from J = 0 to J = 1:")
print(f"   ΔE_shift = -4 * ({D_E:.4e} J) * ({J1+1})³")
print(f"   ΔE_shift = {shift_1_J:.4e} J")
print(f"   In quecto-electronvolts, this is {shift_1_qeV:.4e} qeV\n")

# --- Transition 2: J = 1 to J = 2 ---
J2 = 1
j_term_2 = (J2 + 1)**3
shift_2_J = -4 * D_E * j_term_2
shift_2_qeV = shift_2_J * J_to_eV * eV_to_qeV

print("2. For the transition from J = 1 to J = 2:")
print(f"   ΔE_shift = -4 * ({D_E:.4e} J) * ({J2+1})³")
print(f"   ΔE_shift = {shift_2_J:.4e} J")
print(f"   In quecto-electronvolts, this is {shift_2_qeV:.4e} qeV")