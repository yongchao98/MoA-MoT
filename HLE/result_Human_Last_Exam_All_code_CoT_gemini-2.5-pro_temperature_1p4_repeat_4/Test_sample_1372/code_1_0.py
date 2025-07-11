import math

def solve_star_mass():
    """
    Calculates the mass of a single star in a symmetric triple star system.
    """
    # Given data
    solar_mass_kg = 1.99 * 10**30
    G = 6.67 * 10**-11
    side_length_m = 1.2 * 10**10
    velocity_kms = 125

    # Convert velocity to m/s
    velocity_ms = velocity_kms * 1000

    # The relationship derived from equating gravitational and centripetal forces is:
    # m = v^2 * L / G
    # where m is the mass of one component, v is its tangential velocity,
    # L is the side length of the equilateral triangle, and G is the gravitational constant.

    # Calculate the mass in kilograms
    mass_kg = (velocity_ms**2 * side_length_m) / G

    # Convert the mass to solar masses
    mass_solar = mass_kg / solar_mass_kg

    # Output the explanation and the equation with values
    print("The mass 'm' of a single component is derived from the formula m = v^2 * L / G.")
    print("This comes from equating the net gravitational force F_grav = (G * m^2 / L^2) * sqrt(3) with the centripetal force F_c = m * v^2 / r, where r = L / sqrt(3).")
    print("\nCalculation:")
    print(f"m = ({velocity_ms:.3e})^2 * ({side_length_m:.3e}) / ({G:.3e})")
    print(f"m = {mass_kg:.3e} kg")
    print(f"m = {mass_solar:.3f} solar masses")

    # Final answer rounded to one decimal place
    final_answer = round(mass_solar, 1)
    print(f"\nThe mass of a single component, accurate to one decimal place, is {final_answer} solar masses.")
    print(f"<<<{final_answer}>>>")

solve_star_mass()