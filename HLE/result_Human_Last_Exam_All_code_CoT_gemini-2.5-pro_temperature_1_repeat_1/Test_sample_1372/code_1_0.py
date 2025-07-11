import math

def calculate_star_mass():
    """
    Calculates the mass of a single component in a symmetric triple star system.
    """
    # Given data and constants
    G = 6.67 * 10**-11      # Gravitational constant in m^3 kg^-1 s^-2
    M_solar = 1.99 * 10**30 # Mass of the sun in kg
    L = 1.2 * 10**10        # Side of the equilateral triangle in meters
    T_days = 1.9              # Period in days

    # Step 1: Convert period from days to seconds
    T_sec = T_days * 24 * 60 * 60

    # Step 2: Use the derived formula to calculate the mass in kg.
    # The formula m = (4 * pi^2 * L^3) / (3 * G * T^2) is derived by
    # equating the net gravitational force with the centripetal force.
    
    # Calculate intermediate values to display in the final equation
    pi_val = math.pi
    L_cubed = L**3
    T_sec_sq = T_sec**2
    
    numerator = 4 * (pi_val**2) * L_cubed
    denominator = 3 * G * T_sec_sq
    
    mass_kg = numerator / denominator

    # Step 3: Convert the mass to solar masses
    mass_in_solar_masses = mass_kg / M_solar

    # Step 4: Print the process and the final equation with numbers
    print("The mass (m) of one star is calculated using the formula: m = (4 * pi^2 * L^3) / (3 * G * T^2)\n")
    print("Substituting the values:")
    print(f"L (side length) = {L:.2e} m")
    print(f"G (gravitational constant) = {G:.2e} m^3 kg^-1 s^-2")
    print(f"T (period) = {T_days} days = {T_sec:.2e} s\n")
    
    print("The final equation with the numbers plugged in is:")
    print(f"m = (4 * {pi_val**2:.4f} * {L_cubed:.4e}) / (3 * {G:.2e} * {T_sec_sq:.4e})")
    print(f"m = ({numerator:.4e}) / ({denominator:.4e})\n")

    print("Result:")
    print(f"The mass of a single component is {mass_kg:.4e} kg.")
    print(f"In solar masses, the mass is {mass_in_solar_masses:.1f} M_sun.")

    # Final answer in the required format
    final_answer = round(mass_in_solar_masses, 1)
    print(f"\n<<<{final_answer}>>>")

calculate_star_mass()