import math

def calculate_star_mass():
    """
    Calculates the mass of a single star in a symmetric triple star system.
    """
    # --- Given Data ---
    # Side of the equilateral triangle (m)
    s = 1.2e10
    # Tangential velocity (m/s)
    v = 125 * 1000
    # Gravitational constant (kg^-1 m^3 s^-2)
    G = 6.67e-11
    # Solar mass (kg)
    M_solar = 1.99e30

    # Note: The given period of 1.9 days is inconsistent with the given velocity and geometry.
    # The calculation will proceed using the provided tangential velocity.

    # --- Calculation ---
    # The net gravitational force on one star is F_grav = sqrt(3) * G * m^2 / s^2.
    # The centripetal force is F_c = m * v^2 / r, where r = s / sqrt(3).
    # Equating F_grav and F_c and simplifying leads to: m = v^2 * s / G

    # Calculate the mass in kilograms
    mass_kg = (v**2 * s) / G

    # Convert the mass to solar masses
    mass_solar = mass_kg / M_solar

    # --- Output ---
    print("Equation to find the mass (m) of a single star component:")
    print("m = v^2 * s / G")
    print("\nSubstituting the values:")
    # The instruction asks to output each number in the final equation.
    print(f"m = ({v} m/s)^2 * ({s} m) / ({G} kg^-1 m^3 s^-2)")
    print(f"\nCalculated mass: {mass_kg:.4e} kg")
    print(f"Mass in solar masses: {mass_solar:.4f} M_solar")

    # Final answer rounded to one decimal place
    final_answer = round(mass_solar, 1)
    print(f"\nThe mass of a single component is {final_answer} solar masses.")
    
    # Returning the final answer in the required format
    return final_answer

# Execute the function and print the final answer in the required format
final_result = calculate_star_mass()
print(f"\n<<< {final_result} >>>")
