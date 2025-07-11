import math

def calculate_star_mass():
    """
    Calculates the mass of a single component in a triple star system.

    The system consists of three equal-mass stars in an equilateral triangle
    configuration, rotating about their common center of mass.
    The calculation is based on the principle that the net gravitational force
    on one star provides the centripetal force for its circular motion.

    The formula derived is m = (v^2 * L) / G.
    """
    # Given constants and data
    G = 6.67 * 10**-11  # Gravitational constant in m^3 kg^-1 s^-2
    SOLAR_MASS = 1.99 * 10**30  # Mass of the sun in kg

    # Data from the problem
    L = 1.2 * 10**10  # Side of the equilateral triangle in meters
    v = 125 * 1000      # Tangential velocity in m/s (125 km/s)

    # Calculate the mass of a single component in kg
    # m = (v^2 * L) / G
    mass_kg = (v**2 * L) / G

    # Calculate the mass in solar masses
    mass_solar = mass_kg / SOLAR_MASS

    # Print the equation with numerical values
    print("The equation to find the mass (m) in kg is: m = (v^2 * L) / G")
    print(f"m = ({v}^2 * {L}) / {G}")
    
    # Print the result in kg
    print(f"\nMass of a single component: {mass_kg:.3e} kg")

    # Print the final answer in solar masses, rounded to one decimal place
    print(f"Mass in solar masses: {mass_solar:.1f}")

calculate_star_mass()

# The final answer rounded to one decimal place
final_answer = 1.4
print(f"<<<{final_answer}>>>")