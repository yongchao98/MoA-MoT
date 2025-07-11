import math

# Step 1: Define physical constants
# Reduced Planck's constant in J·s
hbar = 1.054571817e-34
# Atomic mass unit in kg
u = 1.66053906660e-27
# Elementary charge in Coulombs (for eV conversion)
e_charge = 1.602176634e-19

# Step 2: Input data from the problem
# Bond length in meters
r0 = 141.4e-12
# Force constant in N/m
k = 400.0
# Atomic mass of Hydrogen in amu
m_H = 1.008
# Atomic mass of Bromine in amu
m_Br = 79.904

# --- Calculation ---

# Step 3: Calculate the reduced mass (mu) in kg
mu = (m_H * m_Br) / (m_H + m_Br) * u

# Step 4: Calculate the moment of inertia (I)
I = mu * r0**2

# Step 5: Calculate the square of the vibrational energy quantum (hbar * omega)^2 in J^2
# omega = sqrt(k/mu), so (hbar * omega)^2 = hbar^2 * k / mu
hbar_omega_sq = (hbar**2 * k) / mu

# Step 6: Calculate the rotational constant (B) in Joules
B_J = (hbar**2) / (2 * I)

# Step 7: Calculate the centrifugal distortion constant (D) in Joules
D_J = (4 * B_J**3) / hbar_omega_sq

# --- Energy Shifts for Transitions ---

# Step 8.1: Calculate energy shift for J=0 -> J=1 in Joules
# The change in the distortion term is -D * [J_f^2(J_f+1)^2 - J_i^2(J_i+1)^2]
# For J=0 -> J=1, this is -D * [1^2(2)^2 - 0] = -4D
delta_E_J_1 = -4 * D_J

# Step 8.2: Calculate energy shift for J=1 -> J=2 in Joules
# For J=1 -> J=2, this is -D * [2^2(3)^2 - 1^2(2)^2] = -D * [36 - 4] = -32D
delta_E_J_2 = -32 * D_J

# Step 9: Convert energy shifts from Joules to quecto-electronvolts (qeV)
# E_eV = E_J / e_charge
# E_qeV = E_eV / 1e-30 = E_J / (e_charge * 1e-30)
delta_E_qeV_1 = delta_E_J_1 / (e_charge * 1e-30)
delta_E_qeV_2 = delta_E_J_2 / (e_charge * 1e-30)

# Step 10: Print the results as requested.

print("1. Energy shift for transition J=0 to J=1:")
print(f"   ΔE = -4 * D")
print(f"   ΔE = -4 * {D_J:.4e} J = {delta_E_J_1:.4e} J")
print(f"   ΔE (in qeV) = {delta_E_qeV_1:.4e} qeV")
print("")

print("2. Energy shift for transition J=1 to J=2:")
print(f"   ΔE = -32 * D")
print(f"   ΔE = -32 * {D_J:.4e} J = {delta_E_J_2:.4e} J")
print(f"   ΔE (in qeV) = {delta_E_qeV_2:.4e} qeV")
