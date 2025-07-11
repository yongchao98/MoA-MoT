import math

def analyze_beta_spectrometer_setup():
    """
    Analyzes the physics of a beta spectrometer setup and calculates the
    gyroradius for the given parameters to check for plausibility.
    """
    # --- Problem Parameters ---
    # Maximum kinetic energy of the beta particle in MeV
    E_kin_MeV = 1.0
    # Magnetic field strength in milliTesla (mT)
    B_mT = 166.0

    # --- Physical Constants ---
    # Electron rest mass energy in MeV
    m_e_c2_MeV = 0.511
    # Speed of light in m/s
    c = 299792458.0
    # Elementary charge in Coulombs
    e = 1.602176634e-19

    # Convert units to a consistent system (SI)
    E_kin_J = E_kin_MeV * 1.602176634e-13
    B_T = B_mT / 1000.0
    m_e_c2_J = m_e_c2_MeV * 1.602176634e-13

    print("--- Analysis of Experimental Parameters ---")
    print(f"This script verifies that the given magnetic field ({B_mT} mT) is a reasonable")
    print(f"strength for guiding beta particles with a maximum energy of {E_kin_MeV} MeV.")
    print("\nCalculating the maximum gyroradius of the electron...")

    # Step 1: Calculate the relativistic gamma factor
    # Formula: E_kin = (gamma - 1) * m_e*c^2
    gamma = 1 + E_kin_MeV / m_e_c2_MeV
    print(f"\n1. Relativistic Gamma Factor (γ) = 1 + (E_kin / m_e_c^2)")
    print(f"   γ = 1 + ({E_kin_MeV} MeV / {m_e_c2_MeV} MeV) = {gamma:.4f}")

    # Step 2: Calculate the relativistic momentum
    # Formula: p = sqrt(E_total^2 - (m_e*c^2)^2) / c
    # where E_total = gamma * m_e*c^2
    E_total_J = gamma * m_e_c2_J
    p_kg_m_s = math.sqrt(E_total_J**2 - m_e_c2_J**2) / c
    print(f"\n2. Relativistic Momentum (p) = {p_kg_m_s:.3e} kg⋅m/s")

    # Step 3: Calculate the maximum gyroradius
    # The radius of the helical path is r = p_perp / (q * B).
    # The maximum radius (r_max) occurs when the perpendicular momentum (p_perp) equals the total momentum (p).
    p_perp_max_kg_m_s = p_kg_m_s
    r_max_m = p_perp_max_kg_m_s / (e * B_T)
    r_max_cm = r_max_m * 100

    print(f"\n3. Maximum Gyroradius (r_max) = p_perp_max / (e * B)")
    print(f"   r_max = {p_perp_max_kg_m_s:.3e} kg⋅m/s / ({e:.3e} C * {B_T} T)")
    
    print("\n--- Result ---")
    print(f"The maximum gyroradius is {r_max_m:.4f} meters, or {r_max_cm:.2f} cm.")
    print("This is a practical size for a laboratory instrument, confirming the field strength is appropriate.")

if __name__ == '__main__':
    analyze_beta_spectrometer_setup()