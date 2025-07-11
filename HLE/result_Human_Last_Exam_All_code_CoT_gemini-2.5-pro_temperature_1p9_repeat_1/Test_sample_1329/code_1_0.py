import math

def calculate_gyroradius():
    """
    Calculates the maximum gyroradius for a 1 MeV beta particle in a 166 mT magnetic field.
    This demonstrates the feasibility of the field strength proposed in the options.
    """
    # --- Constants ---
    KE_MeV = 1.0  # Kinetic Energy of the electron in Mega-electron Volts
    m_e_c2_MeV = 0.511  # Electron rest mass energy in Mega-electron Volts
    B_tesla = 166e-3  # Magnetic field strength in Tesla (166 mT)
    e_charge = 1.60217663e-19  # Elementary charge in Coulombs
    c_light = 2.99792458e8  # Speed of light in m/s

    # --- Relativistic Momentum Calculation ---
    # 1. Total energy E = KE + m_e*c^2
    E_MeV = KE_MeV + m_e_c2_MeV

    # 2. Relativistic momentum p from E^2 = (pc)^2 + (m_e*c^2)^2
    pc_MeV = math.sqrt(E_MeV**2 - m_e_c2_MeV**2)

    # 3. Convert momentum from MeV/c to SI units (kg*m/s)
    # pc (in Joules) = pc (in MeV) * 1e6 * e
    pc_joules = pc_MeV * 1e6 * e_charge
    # p = pc / c
    p_si = pc_joules / c_light

    # --- Gyroradius Calculation ---
    # The gyroradius is r = p_perp / (q * B). The maximum radius occurs when p_perp = p.
    # q is the electron charge.
    gyroradius_m = p_si / (e_charge * B_tesla)

    # --- Print Results ---
    print("This script verifies if a 166 mT magnetic field can effectively guide a 1 MeV electron.")
    print("We calculate the maximum radius of the electron's spiral path (gyroradius).\n")
    print("--- Calculation Steps ---")
    print(f"1. Total Energy (E) = Kinetic Energy (KE) + Rest Mass Energy (m_e*c^2)")
    print(f"   E = {KE_MeV:.3f} MeV + {m_e_c2_MeV:.3f} MeV = {E_MeV:.3f} MeV\n")

    print(f"2. Relativistic Momentum (p*c) = sqrt(E^2 - (m_e*c^2)^2)")
    print(f"   p*c = sqrt({E_MeV:.3f}^2 - {m_e_c2_MeV:.3f}^2) = {pc_MeV:.3f} MeV\n")

    print(f"3. Gyroradius (r) = p / (q * B)")
    print("   Here, p is the momentum in SI units, q is the electron charge, and B is the magnetic field.")
    print(f"   r = {p_si:.3e} kg*m/s / ({e_charge:.3e} C * {B_tesla:.3f} T)")
    print(f"   r = {gyroradius_m:.4f} m\n")

    print("--- Final Result ---")
    print(f"The maximum gyroradius is {gyroradius_m * 100:.2f} cm.")
    print("This small radius confirms the electron is well-confined, making the proposed field strength effective.")

calculate_gyroradius()