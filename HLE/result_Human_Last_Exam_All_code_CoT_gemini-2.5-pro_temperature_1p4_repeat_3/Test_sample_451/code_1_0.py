import math

# --- Constants ---
# Given values
r0_pm = 141.4      # Bond length in picometers
k_Nm = 400         # Force constant in N/m
m_H_amu = 1.008    # Mass of Hydrogen in amu
m_Br_amu = 79.904  # Mass of Bromine in amu

# Physical constants
h = 6.62607015e-34     # Planck's constant (J·s)
amu_to_kg = 1.660539e-27  # Atomic mass unit to kg conversion
e_charge = 1.602176634e-19 # Electron charge (C)

# Unit conversions
pm_to_m = 1e-12
J_to_eV = 1 / e_charge
eV_to_qeV = 1e30

# --- Calculations ---
# 1. Convert initial values to SI units
r0 = r0_pm * pm_to_m
m_H = m_H_amu * amu_to_kg
m_Br = m_Br_amu * amu_to_kg

# 2. Calculate reduced mass (μ)
mu = (m_H * m_Br) / (m_H + m_Br)

# 3. Calculate moment of inertia (I)
I = mu * r0**2

# 4. Calculate rotational constant (B) in Joules
B = h**2 / (8 * math.pi**2 * I)

# 5. Calculate vibrational frequency (ν_vib) and energy (E_vib)
nu_vib = (1 / (2 * math.pi)) * math.sqrt(k_Nm / mu)
E_vib = h * nu_vib

# 6. Calculate centrifugal distortion constant (D) in Joules
# Formula: D = 4 * B^3 / E_vib^2
D_J = (4 * B**3) / (E_vib**2)

# 7. Calculate energy shifts for transitions in Joules
# Transition 1: J=0 -> J=1; Multiplier is -[1^2*2^2 - 0] = -4
delta_E1_J = -4 * D_J

# Transition 2: J=1 -> J=2; Multiplier is -[2^2*3^2 - 1^2*2^2] = -[36 - 4] = -32
delta_E2_J = -32 * D_J

# 8. Convert energy shifts to quecto-electronvolts (qeV)
J_to_qeV = J_to_eV * eV_to_qeV
delta_E1_qeV = delta_E1_J * J_to_qeV
delta_E2_qeV = delta_E2_J * J_to_qeV

# --- Output Results ---
print("This script calculates the energy shifts due to centrifugal distortion for H-Br.")
print(f"The calculated centrifugal distortion constant is D = {D_J:.4e} J.")
print("-" * 50)

# Output for Transition 1 (J=0 -> J=1)
print("1. Energy shift for rotational transition J = 0 -> J = 1:")
print(f"   Equation: ΔE = -D * [1^2 * (1+1)^2 - 0^2 * (0+1)^2] = -4 * D")
print(f"   Calculation: ΔE = -4 * {D_J:.4e} J = {delta_E1_J:.4e} J")
print(f"   Final Answer: ΔE = {delta_E1_qeV:.3e} qeV")
print("-" * 50)

# Output for Transition 2 (J=1 -> J=2)
print("2. Energy shift for rotational transition J = 1 -> J = 2:")
print(f"   Equation: ΔE = -D * [2^2 * (2+1)^2 - 1^2 * (1+1)^2] = -32 * D")
print(f"   Calculation: ΔE = -32 * {D_J:.4e} J = {delta_E2_J:.4e} J")
print(f"   Final Answer: ΔE = {delta_E2_qeV:.3e} qeV")
