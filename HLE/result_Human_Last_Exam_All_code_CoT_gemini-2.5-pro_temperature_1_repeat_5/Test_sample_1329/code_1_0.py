import math

def calculate_gyroradius():
    """
    This script calculates the maximum gyroradius of a 1 MeV beta particle (electron)
    in a 166 mT magnetic field. This calculation helps to verify that the proposed
    field strength is physically reasonable for guiding electrons in a lab-scale experiment.
    """
    # --- Input Parameters and Constants ---
    T_max_MeV = 1.0  # Maximum kinetic energy in MeV
    B_mT = 166.0     # Magnetic field strength in mT

    # --- Physical Constants ---
    # Electron rest mass energy in MeV
    m_e_c2_MeV = 0.511
    # Elementary charge in Coulombs
    q = 1.60218e-19
    # Speed of light in m/s
    c = 299792458

    # --- Calculations ---
    # Convert units to SI
    T_max_J = T_max_MeV * 1.60218e-13
    m_e_c2_J = m_e_c2_MeV * 1.60218e-13
    B_T = B_mT / 1000.0

    # 1. Calculate total relativistic energy E = T + m_e*c^2
    E_J = T_max_J + m_e_c2_J

    # 2. Calculate relativistic momentum p from E^2 = (pc)^2 + (m_e*c^2)^2
    pc_J = math.sqrt(E_J**2 - m_e_c2_J**2)
    p_kg_m_s = pc_J / c

    # 3. Calculate gyroradius r = p_perp / (q * B)
    # We use the total momentum for p_perp to find the maximum possible radius.
    radius_m = p_kg_m_s / (q * B_T)

    # --- Output Results ---
    print("Verifying the magnetic field strength for guiding a 1 MeV electron.")
    print(f"Maximum Kinetic Energy (T): {T_max_MeV} MeV")
    print(f"Magnetic Field (B): {B_mT} mT")
    print("-" * 30)
    print("The formula for gyroradius is: r = p / (q * B)")
    print("For a 1 MeV electron, the values are:")
    # Final equation with numbers
    print(f"r = {p_kg_m_s:.3e} kg*m/s / ({q:.3e} C * {B_T:.3f} T)")
    print(f"\nThe calculated maximum gyroradius is: {radius_m * 100:.2f} cm.")
    print("This small radius confirms that the field is strong enough to effectively guide the electrons to the detector.")

calculate_gyroradius()