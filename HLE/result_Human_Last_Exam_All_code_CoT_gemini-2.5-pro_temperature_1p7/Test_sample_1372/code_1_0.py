import math

def solve_star_mass():
    """
    Calculates the mass of a single star in a triple star system.
    """
    # Given data
    L = 1.2 * 10**10  # Side of the equilateral triangle in meters
    v_km_s = 125  # Tangential velocity in km/s
    M_sun = 1.99 * 10**30  # Solar mass in kg
    G = 6.67 * 10**-11  # Gravitational constant in kg^-1 m^3 s^-2

    # Convert velocity from km/s to m/s
    v_m_s = v_km_s * 1000

    # Explain the physics and the formula
    print("This problem is solved by equating the net gravitational force on a star to the centripetal force required for its orbit.")
    print("For a star in this equilateral triangle configuration, the relationship is derived as:")
    print("F_gravitational = F_centripetal")
    print("sqrt(3) * G * m^2 / L^2 = m * v^2 / r")
    print("With the orbit radius r = L / sqrt(3), this simplifies to:")
    print("G * m / L = v^2")
    print("\nSolving for the mass (m) of a single star gives the formula:")
    print("m = v^2 * L / G")
    print("\n----------------------------------")

    # Perform the calculation
    print("Calculating the mass with the following values:")
    print(f"Velocity (v) = {v_m_s} m/s")
    print(f"Triangle side length (L) = {L:.1e} m")
    print(f"Gravitational constant (G) = {G:.2e} kg^-1 m^3 s^-2")
    
    # The final equation with numbers
    print("\nThe calculation is:")
    print(f"m = ({v_m_s})^2 * {L:.1e} / {G:.2e}")

    # Calculate mass in kg
    mass_kg = (v_m_s**2 * L) / G
    print(f"\nThe mass of a single component in kg is: {mass_kg:.4e} kg")

    # Convert to solar masses
    mass_solar = mass_kg / M_sun
    print(f"To find the mass in solar masses, we divide by the solar mass ({M_sun:.2e} kg).")
    print(f"\nThe mass of a single component in solar masses is: {mass_solar:.1f}")

    return round(mass_solar, 1)

final_answer = solve_star_mass()
print(f"<<<{final_answer}>>>")
