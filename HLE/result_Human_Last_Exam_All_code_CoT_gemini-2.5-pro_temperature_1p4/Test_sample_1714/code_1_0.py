import math

def calculate_dark_matter_percentage():
    """
    Calculates the dark matter percentage of the Pandora galaxy based on provided data.
    This script simulates the calculations that would be performed on the Wuxing computer.
    """

    # --- Part 1: Define Constants and Initial Values ---
    # Velocity (v): 200 km/s -> 2e5 m/s
    v = 2e5

    # Radius (r): 10 kpc -> 10 * 1000 pc -> 10 * 1000 * 3.086e16 m
    r = 10 * 1000 * 3.086e16

    # Gravitational Constant (G) in m^3 kg^-1 s^-2
    G = 6.674e-11

    # Mass of the Sun (M_sun) in kg
    M_sun = 1.989e30
    
    # --- Part 2: Calculate Luminous Mass (M_lum) ---
    # Luminosity is 2e9 * L_sun
    # Mass/Light ratio is 3 * (M_sun / L_sun)
    # M_lum = (2e9 * L_sun) * (3 * M_sun / L_sun) = 6e9 * M_sun
    lum_factor = 6e9
    m_lum = lum_factor * M_sun

    # --- Part 3: Calculate Total Mass (M_tot) ---
    # Formula: M_tot = (v^2 * r) / G
    v_squared = v**2
    m_tot = (v_squared * r) / G

    # --- Part 4: Calculate Dark Matter Mass and Percentage ---
    m_dark = m_tot - m_lum
    percentage = (m_dark / m_tot) * 100
    
    # --- Part 5: Print the equations with calculated values ---
    print(f"Luminous Mass Equation:")
    print(f"M_lum = Luminosity_Ratio * MassToLight_Ratio * M_sun")
    print(f"M_lum = {lum_factor:.1e} * {M_sun:.4e} kg = {m_lum:.4e} kg\n")

    print(f"Total Mass Equation:")
    print(f"M_tot = (v^2 * r) / G")
    print(f"M_tot = ({v:.1e} m/s)^2 * {r:.4e} m / {G:.4e} m^3 kg^-1 s^-2 = {m_tot:.4e} kg\n")

    print(f"Dark Matter Mass Equation:")
    print(f"M_dark = M_tot - M_lum")
    print(f"M_dark = {m_tot:.4e} kg - {m_lum:.4e} kg = {m_dark:.4e} kg\n")

    print(f"Final Dark Matter Percentage Equation:")
    print(f"Percentage = (M_dark / M_tot) * 100")
    print(f"Percentage = ({m_dark:.4e} kg / {m_tot:.4e} kg) * 100 = {percentage:.1f}%")

calculate_dark_matter_percentage()