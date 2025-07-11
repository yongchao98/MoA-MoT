import math

def solve_star_mass():
    """
    Calculates the mass of a single star in a triple star system
    arranged in an equilateral triangle.
    """

    # --- Given Data ---
    # Side length of the equilateral triangle (m)
    L = 1.2 * 10**10
    # Tangential velocity (km/s)
    v_kms = 125
    # Solar mass (kg)
    M_sun = 1.99 * 10**30
    # Gravitational constant (m^3 kg^-1 s^-2)
    G = 6.67 * 10**-11

    # --- Calculation Steps ---
    # Convert velocity from km/s to m/s
    v_ms = v_kms * 1000

    # Calculate v^2
    v_squared = v_ms**2

    # Calculate the numerator of the equation for mass (v^2 * L)
    numerator = v_squared * L

    # Calculate mass 'm' in kg using the formula m = (v^2 * L) / G
    m_kg = numerator / G

    # Convert mass to solar masses
    m_solar = m_kg / M_sun

    # Round the final answer to one decimal place
    final_answer = round(m_solar, 1)

    # --- Output ---
    print("The mass 'm' of a single component is calculated using the formula: m = (v^2 * L) / G\n")

    print("First, we list the values:")
    print(f"v = {v_kms} km/s = {v_ms:.2e} m/s")
    print(f"L = {L:.2e} m")
    print(f"G = {G:.2e} m^3 kg^-1 s^-2")
    print(f"M_sun = {M_sun:.2e} kg\n")

    print("Now, we substitute the values into the formula to find the mass in kg:")
    print(f"m = (({v_ms:.2e})^2 * {L:.2e}) / {G:.2e}")
    print(f"m = ({v_squared:.4e} * {L:.2e}) / {G:.2e}")
    print(f"m = {numerator:.4e} / {G:.2e}")
    print(f"m = {m_kg:.4e} kg\n")

    print("Finally, we convert the mass from kg to solar masses:")
    print(f"m (solar) = {m_kg:.4e} kg / {M_sun:.2e} kg")
    print(f"m (solar) = {m_solar:.4f} solar masses\n")

    print(f"The mass of a single component, rounded to one decimal place, is {final_answer} solar masses.")

    # Return final answer in specified format
    print(f"<<<{final_answer}>>>")

solve_star_mass()