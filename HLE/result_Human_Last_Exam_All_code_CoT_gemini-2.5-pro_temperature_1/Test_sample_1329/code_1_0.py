import math

def analyze_beta_spectrometer():
    """
    Analyzes the physics of a beta spectrometer and calculates the gyroradius
    for a 1 MeV electron in a 166 mT magnetic field to verify the
    plausibility of the experimental parameters.
    """
    # --- Physical constants and given parameters ---
    K_MeV = 1.0  # Maximum kinetic energy in MeV
    B = 0.166  # Magnetic field in Tesla (166 mT)
    m_e = 9.10938356e-31  # Electron rest mass in kg
    c = 299792458  # Speed of light in m/s
    q = 1.60217662e-19  # Elementary charge in Coulombs
    MeV_to_J = 1.60217662e-13 # Conversion factor from MeV to Joules

    # --- Calculations ---
    # 1. Convert kinetic energy to SI units (Joules)
    K_J = K_MeV * MeV_to_J
    # 2. Calculate electron rest energy in Joules
    E_rest_J = m_e * c**2
    # 3. Calculate the total relativistic energy
    E_total_J = K_J + E_rest_J
    # 4. Calculate relativistic momentum using E_total^2 = (p*c)^2 + (m_e*c^2)^2
    pc_squared = E_total_J**2 - E_rest_J**2
    p = math.sqrt(pc_squared) / c
    # 5. Calculate the maximum gyroradius, r = p / (q * B)
    r = p / (q * B)

    # --- Output the results ---
    print("This script verifies that the proposed magnetic field strength is reasonable for the experiment.")
    print("It calculates the maximum gyroradius for a 1 MeV electron in a 166 mT field.")
    print("-" * 60)

    print("Step 1: Calculate the relativistic momentum (p) of the 1 MeV electron.")
    print(f"The equation for momentum is: p = sqrt(E_total^2 - (m_e*c^2)^2) / c")
    print(f"p = sqrt(({E_total_J:.4e} J)^2 - ({E_rest_J:.4e} J)^2) / {c:.4e} m/s")
    print(f"Calculated momentum p = {p:.4e} kg*m/s")
    print("-" * 60)

    print("Step 2: Calculate the gyroradius (r) in the magnetic field.")
    print("The equation for gyroradius is: r = p / (q * B)")
    print(f"r = {p:.4e} kg*m/s / ({q:.4e} C * {B} T)")
    print(f"Calculated gyroradius r = {r:.4f} meters, or {r*100:.2f} cm")
    print("-" * 60)
    
    print("Conclusion: A maximum gyroradius of a few centimeters is a practical size for a lab experiment,")
    print("confirming the 166 mT field is a sensible value for this application.")

if __name__ == '__main__':
    analyze_beta_spectrometer()
