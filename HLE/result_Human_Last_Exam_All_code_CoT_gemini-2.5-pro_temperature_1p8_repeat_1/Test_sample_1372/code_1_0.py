import math

def calculate_star_mass():
    """
    Calculates the mass of a single component in a symmetric triple star system.
    """
    # Given data and constants
    a = 1.2e10  # Side of the equilateral triangle in meters
    v_kms = 125   # Tangential velocity in km/s
    G = 6.67e-11  # Gravitational constant in m^3 kg^-1 s^-2
    M_sun = 1.99e30  # Solar mass in kg

    # Convert velocity from km/s to m/s
    v = v_kms * 1000

    # The derived formula for the mass of a single star is m = (v^2 * a) / G.
    # We will now plug in the numbers to solve for the mass in kilograms.
    
    # Calculate mass in kilograms
    m_kg = (v**2 * a) / G

    # Convert mass from kilograms to solar masses
    m_solar = m_kg / M_sun

    # Output the steps of the calculation as requested
    print("Step 1: The equation to find the mass of a single star component is derived by equating the net gravitational force with the centripetal force.")
    print("Derived formula: m = (v^2 * a) / G")
    print("\nStep 2: Substitute the given values into the formula.")
    print(f"v (velocity) = {v} m/s")
    print(f"a (triangle side) = {a:.1e} m")
    print(f"G (gravitational constant) = {G:.2e} m^3 kg^-1 s^-2")
    
    print("\nFinal Equation (for mass in kg):")
    # To meet the requirement "output each number in the final equation"
    print(f"m_kg = ({v:.0f}^2 * {a:.1e}) / {G:.2e}")
    print(f"m_kg = {m_kg:.4e} kg")

    print("\nStep 3: Convert the mass from kg to solar masses.")
    print(f"M_sun (solar mass) = {M_sun:.2e} kg")
    print("Final Equation (for mass in solar masses):")
    print(f"m_solar = ({m_kg:.4e}) / ({M_sun:.2e})")
    
    print("\nFinal Answer:")
    # Round the final answer to one decimal place
    print(f"The mass of a single component is {m_solar:.1f} solar masses.")

calculate_star_mass()
<<<1.4>>>