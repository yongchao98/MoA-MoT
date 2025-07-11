import math

# --- Given Constants ---
m_H_amu = 1.008  # Atomic mass of Hydrogen in amu
m_Br_amu = 79.904  # Atomic mass of Bromine in amu
r0_pm = 141.4  # Bond length in picometers
k_Nm = 400.0  # Force constant in N/m

# --- Physical Constants ---
amu_to_kg = 1.660539e-27  # Conversion factor from amu to kg
pm_to_m = 1.0e-12  # Conversion factor from pm to m
h_bar_Js = 1.0545718e-34  # Reduced Planck constant in J*s
e_charge_C = 1.6021766e-19  # Elementary charge in Coulombs

# Step 1: Calculate fundamental molecular properties
# Convert masses to kg
m_H_kg = m_H_amu * amu_to_kg
m_Br_kg = m_Br_amu * amu_to_kg

# Calculate reduced mass (mu) in kg
mu_kg = (m_H_kg * m_Br_kg) / (m_H_kg + m_Br_kg)

# Convert bond length to meters
r0_m = r0_pm * pm_to_m

# Calculate moment of inertia (I) in kg*m^2
I_kg_m2 = mu_kg * r0_m**2

# Step 2: Determine energy constants in Joules
# Calculate rotational constant (B) in Joules
B_J = h_bar_Js**2 / (2 * I_kg_m2)

# Calculate fundamental vibrational energy (hbar*omega) in Joules
omega_energy_J = h_bar_Js * math.sqrt(k_Nm / mu_kg)

# Step 3: Calculate the centrifugal distortion constant (D) in Joules
D_J = (4 * B_J**3) / (omega_energy_J**2)

# Step 6: Convert units and compute final answers
# Conversion factor from Joules to qeV
# 1 J = 1/e_charge_C eV = (1/e_charge_C) * 1e30 qeV
J_to_qeV = (1.0 / e_charge_C) * 1e30

# Calculate D in qeV
D_qeV = D_J * J_to_qeV

# --- Transition 1: J = 0 -> J = 1 ---
J1 = 0
shift1_factor = -4 * (J1 + 1)**3
shift1_qeV = shift1_factor * D_qeV

print("1. Energy shift for the J = 0 -> 1 transition:")
# Output the final equation with all numbers, as requested
print(f"ΔE = -4 * D * (J+1)³ = {shift1_factor} * {D_qeV:.4e} qeV = {shift1_qeV:.4e} qeV")
print("-" * 50)

# --- Transition 2: J = 1 -> J = 2 ---
J2 = 1
shift2_factor = -4 * (J2 + 1)**3
shift2_qeV = shift2_factor * D_qeV

print("2. Energy shift for the J = 1 -> 2 transition:")
# Output the final equation with all numbers, as requested
print(f"ΔE = -4 * D * (J+1)³ = {shift2_factor} * {D_qeV:.4e} qeV = {shift2_qeV:.4e} qeV")
