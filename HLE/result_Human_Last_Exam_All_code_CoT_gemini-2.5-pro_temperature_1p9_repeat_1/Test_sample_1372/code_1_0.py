import math

def solve_star_mass():
    """
    Calculates the mass of a single star in a symmetric triple-star system.

    This function follows a plan to:
    1. Define physical constants and given parameters from the problem.
    2. Note and resolve an inconsistency in the provided data by choosing to use the tangential velocity.
    3. Derive the formula for mass by equating gravitational and centripetal forces.
    4. Substitute the given values into the formula to calculate the mass in kilograms.
    5. Convert the result to solar masses.
    6. Print a step-by-step explanation of the calculation and the final answer rounded to one decimal place.
    """

    # --- Given data and constants ---
    G = 6.67e-11       # Gravitational constant in kg^-1 m^3 s^-2
    M_SOLAR = 1.99e30  # Solar mass in kg
    L = 1.2e10         # Side of the equilateral triangle in m
    V_KMS = 125        # Tangential velocity in km/s
    
    # Convert velocity to m/s
    v_ms = V_KMS * 1000

    print("--- Solving for the Mass of a Single Component Star ---")

    # --- Explanation of inconsistency ---
    print("\nNote on Provided Data:")
    print("The problem provides both a tangential velocity (125 km/s) and a period (1.9 days). These values are inconsistent for the given system geometry.")
    print("Calculating the velocity implied by the period gives a different value (~265 km/s). This solution will proceed using the given tangential velocity of 125 km/s.")
    
    # --- Physics Derivation ---
    print("\nStep 1: Equating Forces")
    print("The net gravitational force on one star from the other two provides the centripetal force that keeps it in orbit.")
    print("F_net_gravity = F_centripetal")
    
    print("\nStep 2: Deriving the Expression for Mass (m)")
    print(" - Net gravitational force: F_net_gravity = (sqrt(3) * G * m^2) / L^2")
    print(" - Orbital radius: r = L / sqrt(3)")
    print(" - Centripetal force: F_centripetal = m * v^2 / r = (m * v^2 * sqrt(3)) / L")
    print("By equating the two force expressions (F_net_gravity = F_centripetal), we can solve for mass m:")
    print("  (sqrt(3) * G * m^2) / L^2 = (m * v^2 * sqrt(3)) / L")
    print("After simplifying, we get a direct relationship:")
    print("  m = (v^2 * L) / G")

    # --- Calculation ---
    mass_kg = (v_ms**2 * L) / G
    mass_solar = mass_kg / M_SOLAR

    print("\nStep 3: Calculating the Mass")
    print("Plugging in the numbers for the final equation m = (v^2 * L) / G:")
    
    # As requested, showing each number in the final equation
    print(f"  m = (({V_KMS} * 1000 m/s)^2 * {L:.1e} m) / {G:.2e} kg^-1 m^3 s^-2")
    print(f"  m = (({v_ms:.2e} m/s)^2 * {L:.1e} m) / {G:.2e} kg^-1 m^3 s^-2")
    print(f"  m = {mass_kg:.4e} kg")

    print("\nStep 4: Converting to Solar Masses")
    print(f"  Mass in solar units = Mass (kg) / Solar Mass (kg)")
    print(f"  Mass = {mass_kg:.4e} kg / {M_SOLAR:.2e} kg")
    print(f"  Mass = {mass_solar:.4f} Solar Masses")

    print("\n--- Final Answer ---")
    print(f"The mass of a single component, accurate to one decimal place, is {mass_solar:.1f} solar masses.")

solve_star_mass()