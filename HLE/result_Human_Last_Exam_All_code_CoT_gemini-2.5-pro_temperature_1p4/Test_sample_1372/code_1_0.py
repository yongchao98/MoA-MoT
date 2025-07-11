import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in an equilateral triple star system.
    """
    # Given data
    side_length_m = 1.2 * 10**10  # Side of the equilateral triangle in meters
    velocity_kms = 125             # Tangential velocity in km/s
    G = 6.67 * 10**-11             # Gravitational constant in kg^-1 m^3 s^-2
    solar_mass_kg = 1.99 * 10**30  # Mass of the sun in kg

    # Convert velocity to m/s
    velocity_ms = velocity_kms * 1000

    # Calculate the mass of a single component in kg
    # Formula: m = v^2 * a / G
    mass_kg = (velocity_ms**2 * side_length_m) / G

    # Convert the mass to solar masses
    mass_solar = mass_kg / solar_mass_kg

    # Output the explanation and the final equation with numbers
    print("The mass 'm' of a single star is determined by equating the net gravitational force to the centripetal force.")
    print("The derived formula is: m = (v^2 * a) / G")
    print("\nCalculation:")
    print(f"mass_kg = ({velocity_ms} m/s)^2 * {side_length_m} m / {G:.2e} kg^-1 m^3 s^-2")
    print(f"mass_kg = {mass_kg:.3e} kg")
    print("\nConverting to solar masses:")
    print(f"mass_solar = {mass_kg:.3e} kg / {solar_mass_kg:.2e} kg")
    print(f"mass_solar = {mass_solar:.2f} Solar Masses")
    
    # Final answer rounded to one decimal place
    final_answer = round(mass_solar, 1)
    
    print("\nThe mass of a single component is:")
    print(f"{final_answer} solar masses.")
    print(f"\n<<<{final_answer}>>>")

calculate_star_mass()
