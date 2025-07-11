import math

# --- Constants ---
KE_eV = 1.0e6  # Kinetic Energy in eV
B = 166e-3     # Magnetic field in Tesla (166 mT)
m_e = 9.10938356e-31  # Electron rest mass in kg
c = 299792458       # Speed of light in m/s
e = 1.60217662e-19  # Elementary charge in Coulombs

# --- Calculations ---

# 1. Convert energies to Joules
KE_J = KE_eV * e
print(f"The kinetic energy of the beta particle is {KE_eV:.1e} eV, which is {KE_J:.4e} J.")

m0c2_J = m_e * c**2
print(f"The rest mass energy of the electron is {m0c2_J:.4e} J.")

# 2. Calculate total relativistic energy
E_total_J = KE_J + m0c2_J
print(f"The total relativistic energy is E = KE + m0c^2 = {E_total_J:.4e} J.")

# 3. Calculate relativistic momentum from the energy-momentum relation: E^2 = (pc)^2 + (m0c^2)^2
# (pc)^2 = E^2 - (m0c^2)^2
pc_squared = E_total_J**2 - m0c2_J**2
p_momentum = math.sqrt(pc_squared) / c
print(f"The relativistic momentum is p = sqrt(E^2 - (m0c^2)^2) / c = {p_momentum:.4e} kg*m/s.")

# 4. Calculate the gyroradius: r = p / (qB)
# This is the radius for a particle with velocity purely perpendicular to B.
# It represents the maximum possible radius in this field.
gyroradius_m = p_momentum / (e * B)
gyroradius_cm = gyroradius_m * 100
print(f"\nFor a magnetic field of {B*1000} mT:")
print(f"The maximum gyroradius is r = p / (qB) = {gyroradius_m:.4f} m, or {gyroradius_cm:.2f} cm.")
print("\nThis radius is small enough to be easily contained in a typical vacuum chamber,")
print("confirming that the proposed field strength is appropriate for guiding the particles.")
