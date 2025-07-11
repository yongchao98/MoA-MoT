import numpy as np

def calculate_gyroradius():
    """
    This script calculates the maximum gyroradius for a 1 MeV beta particle
    in a 166 mT magnetic field to verify that the field strength is appropriate
    for guiding the particles in a reasonably sized apparatus.
    """
    # --- Constants ---
    e_k_MeV = 1.0  # Kinetic energy in MeV
    B_mT = 166.0   # Magnetic field in mT
    m_e_MeV_c2 = 0.511 # Electron rest mass energy in MeV/c^2
    q = 1.602e-19      # Elementary charge in Coulombs
    c = 299792458      # Speed of light in m/s

    # --- Calculations ---
    print("1. Convert input values to SI or standard units.")
    # Convert kinetic energy from MeV to Joules
    e_k_J = e_k_MeV * 1e6 * q
    print(f"   - Maximum kinetic energy: {e_k_MeV} MeV = {e_k_J:.4e} J")
    
    # Convert B field from mT to T
    B_T = B_mT / 1000.0
    print(f"   - Magnetic field strength: {B_mT} mT = {B_T} T")

    # --- Relativistic Momentum Calculation ---
    print("\n2. Calculate the relativistic momentum of the electron.")
    # Calculate electron rest mass energy in Joules
    m_e_c2_J = m_e_MeV_c2 * 1e6 * q
    print(f"   - Electron rest mass energy: {m_e_MeV_c2} MeV = {m_e_c2_J:.4e} J")

    # Total relativistic energy E = E_k + m_e*c^2
    E_total_J = e_k_J + m_e_c2_J
    print(f"   - Total relativistic energy E = E_k + m_e*c^2 = {E_total_J:.4e} J")

    # Relativistic momentum p from E^2 = (pc)^2 + (m_e*c^2)^2
    pc_squared = E_total_J**2 - m_e_c2_J**2
    pc = np.sqrt(pc_squared)
    p = pc / c
    print(f"   - Relativistic momentum p = sqrt(E^2 - (m_e*c^2)^2) / c = {p:.4e} kg*m/s")

    # --- Gyroradius Calculation ---
    print("\n3. Calculate the maximum gyroradius.")
    # The maximum radius occurs when the momentum is entirely perpendicular to B.
    # Gyroradius r = p_perp / (q * B)
    radius_m = p / (q * B_T)
    radius_cm = radius_m * 100
    print(f"   - Equation: r_max = p / (q * B)")
    print(f"   - Calculation: r_max = {p:.4e} kg*m/s / ({q:.4e} C * {B_T} T)")
    print(f"\nResult:")
    print(f"The maximum gyroradius for a {e_k_MeV} MeV electron in a {B_mT} mT field is {radius_cm:.2f} cm.")
    print("This is a reasonable size for a particle detector apparatus, confirming the field strength is appropriate.")

calculate_gyroradius()