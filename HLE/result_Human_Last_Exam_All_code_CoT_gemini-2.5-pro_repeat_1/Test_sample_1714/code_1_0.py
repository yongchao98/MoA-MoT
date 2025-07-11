import math

def calculate_dark_matter_percentage():
    """
    Calculates the percentage of dark matter in the Pandora galaxy,
    using approximations suitable for the Wuxing architecture.
    """

    # I. Define constants based on the problem statement and physical values.
    # The values are approximated to fit into the 'frac' type's limited-size
    # numerator/denominator as described for the Wuxing architecture.
    # For example, G = 6.674e-11 is approximated as 6.7e-11.

    # Gravitational constant G = 6.7e-11 m^3 kg^-1 s^-2
    G = 6.7 * 1e-11

    # Radius r = 10 kpc = 3.086e20 m, approximated as 3.1e20 m
    r = 3.1 * 1e20

    # Velocity v = 200 km/s = 2e5 m/s
    v = 2.0 * 1e5

    # Solar mass M_sun = 1.989e30 kg, approximated as 2.0e30 kg
    m_sun = 2.0 * 1e30

    # Mass-to-Light ratio factor relative to the Sun = 3
    ml_factor = 3.0

    # Luminosity factor relative to the Sun = 2e9
    lum_factor = 2.0 * 1e9

    # II. Calculate Luminous Mass (M_lum)
    # M_lum = (Galaxy's M/L Ratio) * Luminosity
    # M_lum = (3 * M_sun/L_sun) * (2e9 * L_sun) = 3 * 2e9 * M_sun
    M_lum = ml_factor * lum_factor * m_sun

    # III. Calculate Total Mass (M_total) from the velocity curve
    # M_total = (v^2 * r) / G
    M_total = (v**2 * r) / G

    # IV. Calculate the Dark Matter percentage
    # Percentage = ((M_total - M_lum) / M_total) * 100
    percentage = (1 - (M_lum / M_total)) * 100

    # V. Print the final equation and results as requested.
    print(f"Based on the provided data and Wuxing architecture constraints:")
    print(f"Calculated Luminous Mass (M_lum) = {M_lum:.4e} kg")
    print(f"Calculated Total Mass (M_total)   = {M_total:.4e} kg")
    print("\nFinal Equation:")
    print(f"(1 - {M_lum:.4e} / {M_total:.4e}) * 100 = {percentage:.1f}%")

calculate_dark_matter_percentage()