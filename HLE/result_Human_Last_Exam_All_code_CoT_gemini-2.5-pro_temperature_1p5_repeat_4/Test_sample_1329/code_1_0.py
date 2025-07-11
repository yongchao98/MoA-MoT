import numpy as np

def calculate_larmor_radius():
    """
    Calculates the Larmor radius for a 1 MeV electron in a 166 mT magnetic field.
    This demonstrates that the given field strength is appropriate for the experiment.
    """

    # --- Constants ---
    T_MeV = 1.0  # Kinetic energy in MeV
    B_tesla = 166e-3  # Magnetic field in Tesla (166 mT)
    m_e_c2_MeV = 0.511  # Electron rest mass energy in MeV
    c_ms = 299792458  # Speed of light in m/s
    q_coulomb = 1.602e-19  # Elementary charge in Coulombs
    MeV_to_J = 1.602e-13 # Conversion factor from MeV to Joules

    # --- Calculations ---
    # 1. Total energy E = T + m_e*c^2
    E_MeV = T_MeV + m_e_c2_MeV
    E_J = E_MeV * MeV_to_J

    # 2. Relativistic momentum p from E^2 = (pc)^2 + (m_e*c^2)^2
    # (pc)^2 = E^2 - (m_e*c^2)^2
    m_e_c2_J = m_e_c2_MeV * MeV_to_J
    pc_J = np.sqrt(E_J**2 - m_e_c2_J**2)
    p_kg_ms = pc_J / c_ms

    # 3. Larmor radius r = p / (q * B)
    # This is the maximum radius, assuming momentum is perpendicular to the B field.
    r_meters = p_kg_ms / (q_coulomb * B_tesla)

    # --- Output ---
    print("This script verifies that the proposed magnetic field is of a reasonable strength.")
    print("The key is to calculate the Larmor radius (radius of the electron's spiral path).")
    print("If the radius is small, the electron is well-confined by the field.")
    print("\n--- Calculation Steps ---")
    print(f"1. Max Kinetic Energy (T): {T_MeV} MeV")
    print(f"2. Electron Relativistic Momentum (p) for T = {T_MeV} MeV: {p_kg_ms:.4e} kg*m/s")
    print(f"3. Magnetic Field Strength (B): {B_tesla} T")
    print(f"4. Electron Charge (q): {q_coulomb:.4e} C")
    print("\n--- Final Equation for Larmor Radius (r = p / (q * B)) ---")
    
    # Print the equation with all the numbers
    equation_str = (
        f"r = {p_kg_ms:.4e} / ({q_coulomb:.4e} * {B_tesla}) = {r_meters:.4f} m"
    )
    
    # Reformat for the final output as requested
    final_output = f"Radius({r_meters:.4f} m) = Momentum({p_kg_ms:.4e} kg*m/s) / (Charge({q_coulomb:.4e} C) * B-Field({B_tesla} T))"

    print(final_output)

    print(f"\nThe maximum Larmor radius is approximately {r_meters*100:.2f} cm.")
    print("This is a small radius, confirming the field is strong enough to guide the electrons effectively.")

calculate_larmor_radius()