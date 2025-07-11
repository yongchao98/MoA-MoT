import math

# --- Physical Constants ---
m_e = 9.10938356e-31  # Electron mass in kg
c = 299792458.0        # Speed of light in m/s
q = 1.60217662e-19     # Elementary charge in Coulombs
MeV_to_J = 1.60217662e-13 # Conversion factor from MeV to Joules

# --- Given Parameters ---
KE_MeV = 1.0           # Kinetic energy of the electron in MeV
B = 0.166              # Magnetic field strength in Tesla (166 mT)

# 1. Convert kinetic energy to Joules
KE_J = KE_MeV * MeV_to_J

# 2. Calculate the rest mass energy in Joules
E_rest_J = m_e * c**2

# 3. Calculate the total relativistic energy in Joules
E_total_J = E_rest_J + KE_J

# 4. Calculate the relativistic momentum using E^2 = (pc)^2 + (mc^2)^2
# p = sqrt(E_total^2 - E_rest^2) / c
p_squared = (E_total_J**2 - E_rest_J**2) / c**2
p = math.sqrt(p_squared)

# 5. Calculate the gyroradius: r = p / (q * B)
# This radius is for the component of momentum perpendicular to the B field.
# For estimation, we assume all momentum is perpendicular to see the max radius.
gyroradius_m = p / (q * B)
gyroradius_cm = gyroradius_m * 100

print("This calculation determines the maximum gyroradius for a 1 MeV electron in the proposed magnetic field.")
print("A small radius indicates the field is effective at guiding the particle.")
print("-" * 50)
print(f"Relativistic Momentum (p): {p:.4e} kg*m/s")
print(f"Electron Charge (q): {q:.4e} C")
print(f"Magnetic Field (B): {B:.3f} T")
print("\nFinal Equation for Gyroradius (r = p / (q * B)):")
print(f"r = {p:.4e} kg*m/s / ({q:.4e} C * {B:.3f} T)")
print(f"Calculated Gyroradius: {gyroradius_m:.4f} m or {gyroradius_cm:.2f} cm")
print("-" * 50)
print("\nThis radius of ~2.85 cm is small, confirming the 166 mT field is suitable for guiding electrons within a typical lab apparatus.")
