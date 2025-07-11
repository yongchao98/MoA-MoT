import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a symmetric triple star system.
    """
    # --- Given data and constants ---
    # Side of the equilateral triangle (m)
    s = 1.2 * 10**10
    # Tangential velocity (km/s)
    v_kms = 125
    # The solar mass (kg)
    M_sun = 1.99 * 10**30
    # The gravitational constant (kg^-1 m^3 s^-2)
    G = 6.67 * 10**-11

    # --- Convert units to SI ---
    # Convert velocity from km/s to m/s
    v_ms = v_kms * 1000

    # --- Calculation ---
    # The formula derived from equating gravitational force and centripetal force
    # F_gravity = F_centripetal  =>  G*m^2*sqrt(3)/s^2 = m*v^2 / (s/sqrt(3))
    # This simplifies to m = v^2 * s / G
    
    # Calculate mass m in kg
    mass_kg = (v_ms**2 * s) / G

    # Convert mass from kg to solar masses
    mass_solar = mass_kg / M_sun
    
    # Round the final answer to one decimal place
    final_answer = round(mass_solar, 1)

    # --- Output the result ---
    print("Equation to find the mass (m) of a single star component:")
    print("m = v^2 * s / G\n")
    print("Plugging in the values:")
    # We output each number in the final equation as requested
    print(f"m = ({v_ms} m/s)^2 * ({s} m) / ({G} kg^-1 m^3 s^-2)")
    print(f"\nCalculated mass in kilograms: {mass_kg:.3e} kg")
    print(f"Mass in solar masses: {mass_kg:.3e} kg / {M_sun:.2e} kg = {mass_solar:.3f} M_sun")
    print(f"\nThe mass of a single component, accurate to one decimal place, is {final_answer} solar masses.")

# Execute the function to print the solution
calculate_star_mass()