import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a symmetric triple star system.
    """
    # Given data
    G = 6.67e-11  # Gravitational constant in m^3 kg^-1 s^-2
    SOLAR_MASS = 1.99e30  # Solar mass in kg
    L = 1.2e10  # Side of the equilateral triangle in meters
    v_kms = 125  # Tangential velocity in km/s

    # Convert velocity from km/s to m/s
    v_ms = v_kms * 1000

    # Step-by-step derivation explained in the output.
    # The final derived equation for the mass (m) is: m = v^2 * L / G

    # Calculate the mass of a single component in kg
    mass_kg = (v_ms**2 * L) / G

    # Convert the mass to solar masses
    mass_solar = mass_kg / SOLAR_MASS

    # --- Output Section ---
    print("To find the mass of a single star (m), we equate the net gravitational force to the centripetal force.")
    print("The derived equation for the mass is:")
    print("m = (v^2 * L) / G\n")

    print("Using the following values:")
    print(f"  v (tangential velocity) = {v_ms} m/s")
    print(f"  L (triangle side length) = {L:.1e} m")
    print(f"  G (gravitational constant) = {G:.2e} m^3 kg^-1 s^-2\n")

    print("The final equation with the numbers plugged in is:")
    print(f"m = ({v_ms})^2 * {L:.1e} / {G:.2e}")
    
    mass_kg_val = (v_ms)**2 * L / G
    print(f"m = {mass_kg_val:.3e} kg\n")

    print("Converting the mass to solar masses:")
    print(f"Mass in solar masses = {mass_kg_val:.3e} kg / {SOLAR_MASS:.2e} kg")
    
    mass_solar_val = mass_kg_val / SOLAR_MASS
    print(f"Mass in solar masses = {mass_solar_val:.3f}\n")
    
    # Round the final answer to one decimal place
    final_answer = round(mass_solar_val, 1)
    
    print(f"The mass of a single component, accurate to one decimal place, is {final_answer} solar masses.")
    
    # The final answer in the required format
    print(f"<<<{final_answer}>>>")

calculate_star_mass()