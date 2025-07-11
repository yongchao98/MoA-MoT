import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a triple star system.
    """
    # Given data
    L = 1.2e10  # Side length of the equilateral triangle in meters
    v_km_s = 125  # Tangential velocity in km/s
    G = 6.67e-11  # Gravitational constant in kg^-1 m^3 s^-2
    M_sun = 1.99e30  # Solar mass in kg

    # Convert velocity to m/s
    v = v_km_s * 1000

    # The equation to find the mass 'm' of a single component is derived by
    # equating the net gravitational force on one star with the centripetal force.
    # Net gravitational force: F_g = sqrt(3) * G * m^2 / L^2
    # Centripetal force: F_c = m * v^2 / r, where r = L / sqrt(3)
    # Equating F_g = F_c and simplifying gives: m = v^2 * L / G

    print("Equation for the mass of a single component (m):")
    print("m = v^2 * L / G\n")
    print("Substituting the values:")
    print(f"m = ({v} m/s)^2 * ({L:.1e} m) / ({G:.2e} kg^-1 m^3 s^-2)\n")

    # Calculate the mass in kilograms
    mass_kg = (v**2 * L) / G

    print(f"Calculated mass in kilograms: {mass_kg:.4e} kg")

    # Convert mass to solar masses
    mass_solar = mass_kg / M_sun

    print(f"Calculated mass in solar masses: {mass_solar:.4f} M_sun")
    
    # Round to one decimal place for the final answer
    final_answer = round(mass_solar, 1)

    print(f"\nThe mass of a single component, accurate to one decimal place, is {final_answer} solar masses.")
    
    return final_answer

# Run the calculation and store the final answer
final_answer_value = calculate_star_mass()

# Output the final answer in the required format
print(f"\n<<<{final_answer_value}>>>")
