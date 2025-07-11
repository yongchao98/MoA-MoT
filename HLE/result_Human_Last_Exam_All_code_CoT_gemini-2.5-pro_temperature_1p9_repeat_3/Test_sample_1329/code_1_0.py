import math

def calculate_beta_spectrometer_parameters():
    """
    Analyzes the magnetic field requirements for a beta spectrometer.
    This script calculates the gyroradius of a 1 MeV electron in a 166 mT magnetic
    field to demonstrate the feasibility of the setup described in the correct answer.
    """
    # --- Input Parameters & Constants ---
    KE_MeV = 1.0  # Maximum kinetic energy in Mega-electron Volts
    B_T = 166e-3  # Magnetic field strength in Tesla (166 mT)

    # Physical constants
    e = 1.60217663e-19  # Elementary charge in Coulombs
    m0c2_MeV = 0.511  # Electron rest mass energy in MeV
    c = 2.99792458e8  # Speed of light in m/s

    print("To measure the beta spectrum efficiently, we use a magnetic field to guide electrons from the source to the detector.")
    print("The optimal setup uses a gradient, with the field weakest at the source and strongest at the detector.")
    print("Let's calculate the gyroradius (r) of a maximum-energy (1 MeV) electron in the strong field (166 mT) at the detector to ensure it's confined.")
    print("-" * 70)

    # --- Calculation ---
    # 1. Total Energy (E). E = KE + m0c2
    E_MeV = KE_MeV + m0c2_MeV
    E_J = E_MeV * 1e6 * e  # Convert total energy to Joules

    # 2. Relativistic Momentum (p). From E^2 = (pc)^2 + (m0c2)^2
    m0c2_J = m0c2_MeV * 1e6 * e
    pc_J = math.sqrt(E_J**2 - m0c2_J**2)
    p_kgms = pc_J / c

    # 3. Gyroradius (r). From r = p / (|q| * B)
    # For an electron, |q| = e. This is the radius for a particle moving perpendicular to the field.
    gyroradius_m = p_kgms / (e * B_T)
    gyroradius_cm = gyroradius_m * 100

    # --- Output of the Final Equation and its Numbers ---
    print("The final equation to calculate the gyroradius is: r = p / (|q| * B)")
    print("\nWhere the values for a 1 MeV electron are:")
    print(f"  - Relativistic Momentum (p) = {p_kgms:.4e} kg*m/s")
    print(f"  - Electron Charge (|q|)     = {e:.4e} C")
    print(f"  - Magnetic Field (B)        = {B_T} T")

    print(f"\nPlugging in the numbers: r = {p_kgms:.4e} / ({e:.4e} * {B_T})")
    print(f"Resulting Gyroradius (r) = {gyroradius_m:.4f} meters or {gyroradius_cm:.2f} cm.")
    print("-" * 70)
    print("This result shows that the strongest electrons would be confined within a circle of about 2.86 cm radius, which is a practical size for a detector.")

calculate_beta_spectrometer_parameters()