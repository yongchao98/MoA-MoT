import math

def calculate_cable_capacitance():
    """
    Calculates the capacitance of a three-phase cable with a common screen.
    """
    # Step 1: Define given parameters and constants.
    # All lengths are converted from mm to meters for calculation.
    r_w = 11.25 / 1000  # radius of the wires in meters
    m = 17.32 / 1000    # distance from wire center to cable center in meters
    R = 32.32 / 1000    # inner radius of the screen in meters
    epsilon_r = 4.2     # relative permittivity of the insulator
    epsilon_0 = 8.854e-12 # permittivity of free space in F/m

    # Step 2: Calculate the permittivity of the insulator.
    epsilon = epsilon_r * epsilon_0

    # Step 3: Calculate the capacitance in Farads per meter (F/m).
    # The formula for per-phase capacitance is:
    # C = (2 * pi * epsilon) / ln[ (R^3 - m^3) / (3 * r_w * m * R) ]
    try:
        log_term_numerator = R**3 - m**3
        log_term_denominator = 3 * r_w * m * R
        
        if log_term_denominator == 0:
            print("Error: Division by zero in the logarithm term denominator.")
            return
            
        log_argument = log_term_numerator / log_term_denominator
        
        if log_argument <= 0:
            print("Error: The argument for the natural logarithm must be positive.")
            return

        denominator = math.log(log_argument)
        numerator = 2 * math.pi * epsilon
        
        capacitance_F_per_m = numerator / denominator

    except (ValueError, ZeroDivisionError) as e:
        print(f"An error occurred during calculation: {e}")
        return

    # Step 4: Convert the capacitance to microfarads per kilometer (μF/km).
    # Conversion factor from F/m to μF/km is 1e9.
    # (1 F/m = 1e6 μF / 1e-3 km = 1e9 μF/km)
    conversion_factor = 1e9
    capacitance_muF_per_km = capacitance_F_per_m * conversion_factor

    # Step 5: Print the final equation with numbers and the result.
    print("The formula for per-phase capacitance is: C = (2 * π * ε_r * ε_0) / ln[ (R³ - m³) / (3 * r_w * m * R) ]")
    print("\nSubstituting the given values (with lengths in meters):")
    
    # Printing each number in the final equation as requested.
    equation_str = (
        f"C (μF/km) = (2 * {math.pi:.6f} * {epsilon_r} * {epsilon_0}) / "
        f"ln[ ({R}**3 - {m}**3) / (3 * {r_w} * {m} * {R}) ] * {conversion_factor}"
    )
    print(equation_str)

    print(f"\nFinal calculated capacitance: {capacitance_muF_per_km:.4f} μF/km")

    # Returning the final numerical answer for the submission format
    return capacitance_muF_per_km

if __name__ == '__main__':
    result = calculate_cable_capacitance()
    if result is not None:
        print(f"\n<<<C = {result:.4f} μF/km>>>")