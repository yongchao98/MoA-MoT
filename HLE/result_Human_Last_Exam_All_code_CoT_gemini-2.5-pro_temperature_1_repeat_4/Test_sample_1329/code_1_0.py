import math

# --- Constants ---
# Electron rest mass in kg
m_e_kg = 9.10938356e-31
# Speed of light in m/s
c_m_s = 299792458
# Elementary charge in Coulombs
q_C = 1.602176634e-19
# Electron rest energy in MeV
m_e_c2_MeV = 0.511

# --- Given Parameters ---
# Kinetic energy of the beta particle in MeV
T_MeV = 1.0
# Magnetic field strength in Tesla
B_T = 166e-3  # 166 mT

# --- Calculations ---
# 1. Calculate the relativistic gamma factor
# T = (gamma - 1) * m_e * c^2
gamma = 1 + T_MeV / m_e_c2_MeV

# 2. Calculate the particle's velocity
# gamma = 1 / sqrt(1 - v^2/c^2) => v = c * sqrt(1 - 1/gamma^2)
v_fraction_c = math.sqrt(1 - 1 / (gamma**2))
v_m_s = v_fraction_c * c_m_s

# 3. Calculate the relativistic momentum
# p = gamma * m_e * v
p_kg_m_s = gamma * m_e_kg * v_m_s

# 4. Calculate the maximum Larmor radius (for p_perpendicular = p)
# r = p / (q * B)
r_m = p_kg_m_s / (q_C * B_T)
r_cm = r_m * 100

# --- Output the results and the final equation ---
print("This calculation determines the maximum Larmor radius for a 1.0 MeV electron in the given magnetic field.")
print("This helps verify that the proposed field strength is appropriate for a lab-scale experiment.")
print("-" * 60)
print(f"Given values:")
print(f"  Kinetic Energy (T) = {T_MeV} MeV")
print(f"  Magnetic Field (B) = {B_T * 1000} mT")
print("-" * 60)
print("Calculation steps:")
print(f"1. Relativistic Gamma (γ):")
print(f"   γ = 1 + T_MeV / m_e_c2_MeV = 1 + {T_MeV:.1f} / {m_e_c2_MeV:.3f} = {gamma:.4f}")
print()
print(f"2. Velocity (v):")
print(f"   v = c * sqrt(1 - 1/γ^2) = c * sqrt(1 - 1/{gamma:.4f}^2) = {v_fraction_c:.4f}c")
print()
print(f"3. Relativistic Momentum (p):")
print(f"   p = γ * m_e * v = {gamma:.4f} * {m_e_kg:.3e} kg * {v_m_s:.3e} m/s = {p_kg_m_s:.3e} kg*m/s")
print()
print("4. Maximum Larmor Radius (r_max):")
print("   The final equation is r_max = p / (q * B)")
print("   Plugging in the numbers:")
print(f"   r_max = {p_kg_m_s:.3e} kg*m/s / ({q_C:.3e} C * {B_T:.3f} T)")
print(f"\nResult:")
print(f"   The maximum Larmor radius is {r_cm:.2f} cm.")
print("-" * 60)
print("This radius is a manageable size for a laboratory apparatus, confirming")
print("the feasibility of the proposed magnetic field strength.")