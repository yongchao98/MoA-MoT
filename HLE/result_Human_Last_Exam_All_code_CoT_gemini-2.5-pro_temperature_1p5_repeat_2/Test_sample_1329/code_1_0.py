import math

# Description of the calculation
print("This script calculates the cyclotron radius of a 1 MeV beta particle (electron) in a 166 mT magnetic field.")
print("This helps to understand how the magnetic field confines the particles to guide them towards the detector.\n")

# --- Constants ---
# Kinetic energy of the electron in MeV
T_MeV = 1.0
# Magnetic field strength in Tesla
B_T = 166e-3  # 166 mT
# Electron rest mass energy in MeV/c^2
m0c2_MeV = 0.511
# Speed of light in m/s
c_ms = 299792458
# Elementary charge in Coulombs
q_C = 1.602176634e-19
# Joules per MeV
J_per_MeV = 1.602176634e-13

# --- Calculation ---
# 1. Calculate the total energy of the electron
E_MeV = T_MeV + m0c2_MeV
print(f"The kinetic energy T is {T_MeV} MeV.")
print(f"The electron rest mass energy m0*c^2 is {m0c2_MeV} MeV.")
print(f"The total relativistic energy E = T + m0*c^2 is {E_MeV:.3f} MeV.\n")

# 2. Calculate the relativistic momentum of the electron
# From E^2 = (p*c)^2 + (m0*c^2)^2
pc_MeV = math.sqrt(E_MeV**2 - m0c2_MeV**2)
print(f"The momentum in units of MeV/c is p*c = sqrt(E^2 - (m0*c^2)^2) = {pc_MeV:.3f} MeV.")

# 3. Convert momentum to SI units (kg*m/s)
p_si = (pc_MeV * J_per_MeV) / c_ms
print(f"This momentum in SI units is p = {p_si:.3e} kg*m/s.\n")

# 4. Calculate the cyclotron radius
# The formula for the radius of the path is r = p / (q * B)
# This assumes the momentum is perpendicular to the B field (worst-case for confinement)
radius_m = p_si / (q_C * B_T)
radius_cm = radius_m * 100

# --- Output the final result ---
print(f"The magnetic field B is {B_T*1000} mT.")
print(f"The charge of an electron q is {q_C:.3e} C.")
print(f"The resulting maximum cyclotron radius is r = p / (q * B)")
print(f"r = {p_si:.3e} / ({q_C:.3e} * {B_T}) = {radius_m:.4f} meters")
print(f"\nFinal Result: The maximum radius of the electron's spiral path is {radius_cm:.2f} cm.")
print("This small radius shows that the field is strong enough to effectively guide the particles.")