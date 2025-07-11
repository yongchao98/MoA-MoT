import math

def solve_star_mass():
    """
    Calculates the mass of a single star in a triple star system.
    """
    # --- Given data and constants ---

    # Side of the equilateral triangle (m)
    a = 1.2 * 10**10
    # Tangential velocity (km/s), converted to m/s
    v = 125 * 1000
    # Gravitational constant (m^3 kg^-1 s^-2)
    G = 6.67 * 10**-11
    # Solar mass (kg)
    M_solar = 1.99 * 10**30
    
    # --- Explanation and Derivation ---
    
    print("This problem is solved by equating the net gravitational force on one star with the centripetal force required for its circular orbit.")
    print("The given period of 1.9 days is inconsistent with the given velocity and distance, so we proceed using the instantaneous velocity.")
    print("\nThe net gravitational force on one star from the other two is F_g = (G * m^2 / a^2) * sqrt(3).")
    print("The centripetal force is F_c = m * v^2 / r, where the orbit radius r = a / sqrt(3).")
    print("Equating F_g = F_c and simplifying leads to the formula for the mass 'm' in kilograms:")
    print("m = (v^2 * a) / G")
    print("\nTo find the answer in solar masses, we divide by the solar mass M_solar.")
    
    # --- Calculation ---
    
    # The final equation for mass in solar masses is: (v^2 * a) / (G * M_solar)
    print("\nThe final equation with the numerical values is:")
    print(f"Mass in solar masses = ({v:.3e}^2 * {a:.3e}) / ({G:.3e} * {M_solar:.3e})")
    
    # Calculate the numerator and denominator separately for clarity
    numerator = v**2 * a
    denominator = G * M_solar
    
    print(f"Mass in solar masses = {numerator:.3e} / {denominator:.3e}")
    
    # Calculate the final result
    mass_in_solar_masses = numerator / denominator
    
    # Round to one decimal place as requested
    final_answer = round(mass_in_solar_masses, 1)
    
    print(f"\nResult: {mass_in_solar_masses:.4f} solar masses")
    print(f"The mass of a single component, accurate to one decimal place, is {final_answer} solar masses.")
    
    # Return the final answer in the required format
    return final_answer

# Execute the function and capture the result
final_result = solve_star_mass()
print(f"\n<<< {final_result} >>>")
