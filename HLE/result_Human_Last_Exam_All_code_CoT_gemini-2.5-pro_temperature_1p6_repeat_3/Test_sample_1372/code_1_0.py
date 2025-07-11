import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a symmetric triple star system.
    """
    # Step 1: Define physical constants and given parameters in SI units.
    # Side of the equilateral triangle (m)
    a = 1.2e10
    # Tangential velocity (m/s), converted from 125 km/s
    v = 125e3
    # Gravitational constant (m^3 kg^-1 s^-2)
    G = 6.67e-11
    # Solar mass (kg)
    M_sun = 1.99e30

    # Step 2: Calculate the mass 'm' in kilograms using the derived formula.
    # The formula m = (a * v^2) / G is derived by equating the centripetal force
    # with the net gravitational force on one of the stars.
    mass_kg = (a * v**2) / G

    # Step 3: Convert the mass from kilograms to solar masses.
    mass_solar = mass_kg / M_sun

    # Step 4: Print the calculation steps for clarity.
    print("The mass 'm' of a single star is calculated using the formula: m = (a * v^2) / G")
    print("\nSubstituting the given values into the formula:")
    print(f"m = ({a} * {v**2}) / {G}")
    print(f"m = {mass_kg:.4e} kg")

    print("\nConverting this mass to solar masses:")
    print(f"m_solar = mass_kg / M_sun")
    print(f"m_solar = {mass_kg:.4e} kg / {M_sun:.2e} kg")
    
    # Step 5: Present the final answer rounded to one decimal place.
    final_answer = round(mass_solar, 1)
    print(f"\nThe mass of a single component is approximately {final_answer} solar masses.")

calculate_star_mass()