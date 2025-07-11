import math

def calculate_larmor_radius():
    """
    Calculates the Larmor radius for a 1 MeV electron in a 166 mT magnetic field.
    This demonstrates the physical feasibility of the parameters given in the options.
    """
    # Physical constants
    m_e = 9.1093837e-31  # Electron rest mass in kg
    c = 2.99792458e8     # Speed of light in m/s
    q = 1.60217663e-19   # Elementary charge in Coulombs
    
    # Given parameters
    kinetic_energy_MeV = 1.0  # Kinetic energy in MeV
    magnetic_field_T = 0.166  # Magnetic field in Tesla (166 mT)

    # --- Calculations ---
    # Convert energy to Joules
    MeV_to_J = 1.60217663e-13
    kinetic_energy_J = kinetic_energy_MeV * MeV_to_J
    
    # Electron rest mass energy in Joules
    rest_energy_J = m_e * c**2
    
    # Total relativistic energy in Joules
    total_energy_J = rest_energy_J + kinetic_energy_J
    
    # Calculate relativistic momentum from the energy-momentum relation: E^2 = (pc)^2 + (m_0c^2)^2
    # (pc)^2 = E^2 - (m_0c^2)^2
    pc_squared = total_energy_J**2 - rest_energy_J**2
    pc = math.sqrt(pc_squared)
    momentum = pc / c
    
    # Calculate Larmor radius: r = p / (q * B)
    # This is the maximum radius, assuming the electron is emitted perpendicular to the B-field.
    larmor_radius_m = momentum / (q * magnetic_field_T)
    
    print("This calculation checks if the proposed magnetic field is reasonable for a lab experiment.")
    print("We calculate the maximum spiral radius for a 1 MeV electron.\n")
    print("Final Equation: Larmor Radius (r) = Momentum (p) / (Charge (q) * Magnetic Field (B))")
    print("-" * 50)
    print(f"Relativistic momentum (p): {momentum:.4g} kg*m/s")
    print(f"Elementary charge (q): {q:.4g} C")
    print(f"Magnetic field (B): {magnetic_field_T:.4g} T")
    print("-" * 50)
    print(f"Calculated max Larmor radius: {larmor_radius_m:.4f} m")
    print(f"Calculated max Larmor radius: {larmor_radius_m * 100:.2f} cm")
    print("\nA radius of a few centimeters is a manageable size for a particle detector setup, confirming the feasibility of the chosen parameters.")

calculate_larmor_radius()