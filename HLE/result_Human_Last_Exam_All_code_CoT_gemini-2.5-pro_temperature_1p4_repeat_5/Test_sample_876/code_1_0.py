import math

def calculate_cable_capacitance():
    """
    Calculates the capacitance of a three-phase cable with a common screen.
    """
    # Step 1: Define constants and variables in SI units
    r_w = 11.25 / 1000  # radius of the wires in meters
    m = 17.32 / 1000    # distance from wire center to cable center in meters
    R = 32.32 / 1000    # external radius of the cable in meters
    epsilon_r = 4.2     # relative permittivity of the insulator (dimensionless)
    epsilon_0 = 8.854e-12 # permittivity of free space in F/m

    # Combined permittivity
    epsilon = epsilon_r * epsilon_0

    # Step 2: Calculate the terms of the formula
    # Numerator of the main formula
    formula_numerator = 2 * math.pi * epsilon

    # Argument of the natural logarithm
    log_arg_numerator = R**3 - m**3
    log_arg_denominator = 3 * m * r_w * R
    log_arg = log_arg_numerator / log_arg_denominator

    # Denominator of the main formula
    formula_denominator = math.log(log_arg)

    # Step 3: Calculate capacitance in Farads per meter (F/m)
    capacitance_F_per_m = formula_numerator / formula_denominator

    # Step 4: Convert capacitance to microfarads per kilometer (uF/km)
    # 1 F/m = 1e6 uF / 1e-3 km = 1e9 uF/km
    capacitance_muF_per_km = capacitance_F_per_m * 1e9

    # Step 5: Print the detailed calculation process and result
    print("The capacitance per phase (C) for a three-phase cable with a common screen is calculated using the formula:")
    print("C = (2 * pi * ε₀ * εᵣ) / ln((R³ - m³) / (3 * m * r_w * R))")
    print("\nSubstituting the given values into the equation (in SI units):")
    
    # Print the equation with all the numbers
    equation_str = (
        f"C (F/m) = (2 * {math.pi:.5f} * {epsilon_0:.3e} F/m * {epsilon_r}) / "
        f"ln((({R} m)³ - ({m} m)³) / (3 * {m} m * {r_w} m * {R} m))"
    )
    print(equation_str)
    
    # Print the final result
    print(f"\nAfter performing the calculation:")
    print(f"The capacitance per phase of the cable is {capacitance_muF_per_km:.4f} \u03BCF/km.")
    
    # Return the final numerical answer for the platform
    return capacitance_muF_per_km

if __name__ == '__main__':
    final_capacitance = calculate_cable_capacitance()
    # The final answer format for the platform
    # print(f"<<<{final_capacitance:.4f}>>>")
    # For this exercise, the final answer will be embedded in the output.
    pass
    print(f"\n<<<{final_capacitance:.4f}>>>")
