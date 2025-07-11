import math

# Step 1: Define Constants and Parameters
# Physical constants
h = 6.62607015e-34      # Planck's constant (J·s)
amu_to_kg = 1.660539e-27 # Atomic mass unit to kg conversion
e_charge = 1.602176634e-19 # Elementary charge (C)
pi = math.pi

# Input parameters for H-Br molecule
m_H = 1.008             # Mass of Hydrogen (amu)
m_Br = 79.904           # Mass of Bromine (amu)
r0 = 141.4e-12          # Bond length (m)
k = 400                 # Force constant (N/m)

# Step 2: Calculate Reduced Mass (mu) and Moment of Inertia (I)
mu_amu = (m_H * m_Br) / (m_H + m_Br)
mu_kg = mu_amu * amu_to_kg
I = mu_kg * r0**2

# Step 3: Calculate Spectroscopic Constants (B, nu, D)
# Rotational constant in Joules
B_joules = h**2 / (8 * pi**2 * I)
# Vibrational frequency in Hz
nu_hz = (1 / (2 * pi)) * math.sqrt(k / mu_kg)
# Centrifugal distortion constant in Joules
D_joules = (4 * B_joules**3) / (h * nu_hz)**2

print(f"The energy shift due to centrifugal distortion in a transition J -> J+1 is given by ΔE = -4 * D * (J+1)^3.")
print(f"Calculated Centrifugal Distortion Constant, D = {D_joules:.4e} J")
print("-" * 50)

# Step 4 & 5: Calculate, Convert, and Print Energy Shifts

# Define the conversion factor from Joules to quecto-electronvolts (qeV)
qeV_in_eV = 1e-30
J_per_qeV = e_charge * qeV_in_eV
J_to_qeV = 1 / J_per_qeV

# -- Calculation for Transition 1: J = 0 to J = 1 --
J1 = 0
# The equation for the energy shift is ΔE_shift = -4 * D * (J+1)^3
delta_E1_J = -4 * D_joules * (J1 + 1)**3
delta_E1_qeV = delta_E1_J * J_to_qeV

print("1. For the transition from J = 0 to J = 1:")
print(f"   Equation: ΔE_shift = -4 * D * (J+1)^3")
print(f"   Calculation: ΔE_shift = -4 * ({D_joules:.4e}) * ({J1}+1)^3 = {delta_E1_J:.4e} J")
print(f"   The resulting energy shift is: {delta_E1_qeV:.3e} qeV")
print("-" * 50)


# -- Calculation for Transition 2: J = 1 to J = 2 --
J2 = 1
# The equation for the energy shift is ΔE_shift = -4 * D * (J+1)^3
delta_E2_J = -4 * D_joules * (J2 + 1)**3
delta_E2_qeV = delta_E2_J * J_to_qeV

print("2. For the transition from J = 1 to J = 2:")
print(f"   Equation: ΔE_shift = -4 * D * (J+1)^3")
print(f"   Calculation: ΔE_shift = -4 * ({D_joules:.4e}) * ({J2}+1)^3 = {delta_E2_J:.4e} J")
print(f"   The resulting energy shift is: {delta_E2_qeV:.3e} qeV")
print("-" * 50)


final_answer = f"For J=0->1: {delta_E1_qeV:.3e} qeV; For J=1->2: {delta_E2_qeV:.3e} qeV."
print(f"\nFinal Answer Summary:\n{final_answer}")
<<<For J=0->1: -1.762e+23 qeV; For J=1->2: -1.410e+24 qeV.>>>