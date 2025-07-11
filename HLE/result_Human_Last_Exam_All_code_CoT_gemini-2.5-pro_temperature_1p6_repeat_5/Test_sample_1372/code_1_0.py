import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a triple star system.
    """
    # Given data
    L = 1.2 * 10**10  # Side of the equilateral triangle in meters
    v_km_s = 125      # Tangential velocity in km/s
    
    # Constants
    G = 6.67 * 10**-11  # Gravitational constant in kg^-1 m^3 s^-2
    M_sun = 1.99 * 10**30  # Solar mass in kg

    # Convert velocity to m/s
    v = v_km_s * 1000

    # The physics:
    # The net gravitational force on one star from the other two provides the centripetal force.
    # F_net = F_centripetal
    # sqrt(3) * G * M^2 / L^2 = M * v^2 / r
    # The radius of orbit r (distance from vertex to centroid) is r = L / sqrt(3).
    # sqrt(3) * G * M^2 / L^2 = M * v^2 / (L / sqrt(3))
    # After simplification, we get the formula for mass M:
    # M = v^2 * L / G

    # Calculate the mass of a single component in kg
    mass_kg = (v**2 * L) / G

    # Convert the mass to solar masses
    mass_solar = mass_kg / M_sun

    # Print the equation with the numbers used
    print("Derived formula for mass M: M = v^2 * L / G")
    print(f"M = ({v:.0f} m/s)^2 * {L:.1e} m / {G:.2e} kg^-1 m^3 s^-2")
    
    # Print the intermediate and final results
    print(f"Mass of one component: {mass_kg:.4e} kg")
    print(f"Mass in solar masses: {mass_kg:.4e} kg / {M_sun:.2e} kg = {mass_solar:.3f} Solar Masses")
    print(f"Rounded to one decimal place, the mass is: {mass_solar:.1f} Solar Masses")


calculate_star_mass()
<<<1.4>>>