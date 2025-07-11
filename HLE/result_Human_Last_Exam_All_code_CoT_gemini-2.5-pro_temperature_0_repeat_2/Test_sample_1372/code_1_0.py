import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a symmetric triple star system.
    """
    # --- Given Data ---
    # Gravitational constant in kg^-1 m^3 s^-2
    G = 6.67 * 10**-11
    # Side of the equilateral triangle in meters
    L = 1.2 * 10**10
    # Tangential velocity in m/s (125 km/s)
    v = 125 * 1000
    # Solar mass in kg
    M_solar = 1.99 * 10**30

    # --- Derivation Explanation ---
    # The net gravitational force on one star from the other two is F_g = sqrt(3) * G * m^2 / L^2.
    # The radius of the orbit (distance from a vertex to the centroid of the triangle) is r = L / sqrt(3).
    # The centripetal force is F_c = m * v^2 / r = sqrt(3) * m * v^2 / L.
    # By equating the gravitational force and the centripetal force (F_g = F_c),
    # we get: sqrt(3) * G * m^2 / L^2 = sqrt(3) * m * v^2 / L.
    # This simplifies to the final equation for mass (m): m = v^2 * L / G.

    print("The final equation for the mass (m) of a single component is: m = v^2 * L / G")
    print("\nPlugging in the values:")
    # Output each number in the final equation
    print(f"m = ({v:.1f} m/s)^2 * ({L:.1e} m) / ({G:.2e} kg^-1 m^3 s^-2)")

    # --- Calculation ---
    # Calculate the mass of a single star in kilograms
    mass_kg = (v**2 * L) / G

    # Convert the mass to solar masses
    mass_in_solar_masses = mass_kg / M_solar

    print(f"\nThe calculated mass of a single component is {mass_kg:.4e} kg.")
    print(f"This is equivalent to {mass_in_solar_masses:.1f} solar masses.")

    # --- Final Answer ---
    # The final answer is formatted as requested.
    print(f"\n<<<{mass_in_solar_masses:.1f}>>>")

calculate_star_mass()