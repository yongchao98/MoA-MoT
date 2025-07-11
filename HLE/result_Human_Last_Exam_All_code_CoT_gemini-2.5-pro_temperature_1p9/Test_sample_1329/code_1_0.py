import math

def calculate_larmor_radius():
    """
    Calculates the Larmor radius for a 1 MeV electron in a 166 mT magnetic field.
    This demonstrates the physical scale of the particle's spiral path.
    """
    # Constants
    K_MeV = 1.0  # Kinetic Energy in MeV
    B_T = 166e-3  # Magnetic Field in Tesla (166 mT)
    
    m_e_MeV = 0.511  # Electron rest mass in MeV/c^2
    m_e_kg = 9.1093837e-31  # Electron rest mass in kg
    c = 299792458.0  # Speed of light in m/s
    q = 1.60217663e-19  # Elementary charge in Coulombs

    print(f"Calculating Larmor radius for an electron with:")
    print(f"Kinetic Energy (K) = {K_MeV} MeV")
    print(f"Magnetic Field (B) = {B_T} T\n")

    # --- Relativistic Calculations ---
    
    # 1. Total Energy (E = K + E_0)
    E_total_MeV = K_MeV + m_e_MeV
    print(f"Total Relativistic Energy (E) = {K_MeV} MeV + {m_e_MeV} MeV = {E_total_MeV:.3f} MeV")

    # 2. Lorentz Factor (gamma = E / E_0)
    gamma = E_total_MeV / m_e_MeV
    print(f"Lorentz Factor (gamma) = {E_total_MeV:.3f} MeV / {m_e_MeV} MeV = {gamma:.3f}")

    # 3. Velocity (v = c * sqrt(1 - 1/gamma^2))
    v = c * math.sqrt(1 - 1/gamma**2)
    print(f"Velocity (v) = {v:.3e} m/s (or {v/c:.3f}c)")

    # 4. Relativistic Momentum (p = gamma * m * v)
    # We calculate momentum for the case where the electron is emitted perpendicular
    # to the B-field, which gives the maximum Larmor radius.
    p_perp = gamma * m_e_kg * v
    print(f"Maximum Perpendicular Momentum (p_perp) = {p_perp:.3e} kg*m/s\n")
    
    # 5. Larmor Radius (r = p_perp / (q * B))
    radius_m = p_perp / (q * B_T)

    print("--- Larmor Radius Calculation ---")
    print("The formula for Larmor radius is: r = p_perp / (q * B)")
    print("Plugging in the numbers:")
    # The prompt requires printing the equation with the final numbers
    print(f"r = {p_perp:.4e} kg*m/s / ({q:.4e} C * {B_T} T)")
    print(f"Maximum Larmor Radius (r) = {radius_m:.4f} meters")
    print(f"Maximum Larmor Radius (r) = {radius_m * 100:.2f} cm")
    
    print("\nThis radius of a few centimeters is a manageable size for a laboratory detector system.")

calculate_larmor_radius()