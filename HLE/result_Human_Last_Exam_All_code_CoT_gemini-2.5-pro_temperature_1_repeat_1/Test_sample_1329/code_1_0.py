import math

def analyze_beta_spectrometry_setup():
    """
    Calculates the Larmor radius for a 1 MeV electron in a 166 mT magnetic field
    to verify the feasibility of the proposed experimental setup.
    """
    # --- Constants and Parameters ---
    # Physical constants
    q = 1.60217663e-19  # Elementary charge in Coulombs
    m_e_c2_MeV = 0.511  # Electron rest mass energy in MeV
    c = 299792458       # Speed of light in m/s
    MeV_to_J = 1.60217663e-13 # Conversion factor from MeV to Joules

    # Experiment parameters from the problem description
    T_max_MeV = 1.0       # Maximum kinetic energy of beta particle in MeV
    B_T = 166e-3          # Magnetic field strength in Tesla (166 mT)

    # --- Explanation and Calculation ---
    print("To assess the proposed magnetic field, we calculate the Larmor radius (gyroradius)")
    print("for the most energetic electrons to ensure they are effectively guided.")
    
    print("\nStep 1: Calculate the total energy (E) of the electron.")
    E_MeV = T_max_MeV + m_e_c2_MeV
    print(f"E = Kinetic Energy + Rest Mass Energy")
    print(f"E = {T_max_MeV:.3f} MeV + {m_e_c2_MeV:.3f} MeV = {E_MeV:.3f} MeV")

    print("\nStep 2: Calculate the electron's momentum (p) using the relativistic relation E^2 = (p*c)^2 + (m_e*c^2)^2.")
    pc_MeV = math.sqrt(E_MeV**2 - m_e_c2_MeV**2)
    print(f"p*c = sqrt(E^2 - (m_e*c^2)^2)")
    print(f"p*c = sqrt({E_MeV:.3f}^2 - {m_e_c2_MeV:.3f}^2) MeV = {pc_MeV:.3f} MeV")

    print("\nStep 3: Convert momentum from 'MeV/c' to SI units (kg*m/s).")
    p_si = (pc_MeV * MeV_to_J) / c
    print(f"p = (p*c [in J]) / c")
    print(f"p = ({pc_MeV * MeV_to_J:.3e} J) / {c:.3e} m/s = {p_si:.3e} kg*m/s")

    print("\nStep 4: Calculate the Larmor radius (r_L) for the maximum case (momentum perpendicular to B).")
    r_L = p_si / (q * B_T)
    print(f"r_L = p / (|q| * B)")
    print(f"r_L = {p_si:.3e} kg*m/s / ({q:.3e} C * {B_T:.3f} T)")
    
    r_L_cm = r_L * 100
    print(f"\nThe calculated maximum Larmor radius is: {r_L_cm:.2f} cm.")
    print("\nThis radius is small enough for effective particle guiding in a typical lab setup, confirming the field strength is appropriate.")

analyze_beta_spectrometry_setup()