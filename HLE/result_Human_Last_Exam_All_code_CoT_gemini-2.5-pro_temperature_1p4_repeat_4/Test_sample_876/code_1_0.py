import math

def calculate_capacitance():
    """
    Calculates and prints the capacitance of a three-phase cable with a common screen.
    """
    # Step 1: Define parameters and constants
    # Given values
    r_w_mm = 11.25  # radius of wires in mm
    m_mm = 17.32    # distance from wire center to cable center in mm
    R_mm = 32.32    # inner radius of the screen in mm
    epsilon_r = 4.2 # relative permittivity of the insulator

    # Convert mm to m for consistent calculations
    r_w = r_w_mm / 1000.0
    m = m_mm / 1000.0
    R = R_mm / 1000.0

    # Physical constants
    epsilon_0 = 8.854e-12 # Permittivity of free space in F/m

    # Step 2: Calculate absolute permittivity
    epsilon = epsilon_0 * epsilon_r

    # Step 3: Calculate the capacitance using the formula step-by-step
    # This section breaks down the formula for clarity in the output
    
    # Numerator of the main formula
    numerator = 2 * math.pi * epsilon

    # Components of the denominator argument
    m_div_R = m / R
    m_div_R_pow6 = m_div_R ** 6
    term1_in_root = 1 - m_div_R_pow6
    term2_root = term1_in_root ** (1/3)
    R_div_rw = R / r_w
    ln_argument = R_div_rw * term2_root
    
    # Final denominator
    denominator = math.log(ln_argument)

    # Capacitance in Farads per meter (F/m)
    C_F_per_m = numerator / denominator

    # Step 4: Convert the result to microfarads per kilometer (uF/km)
    C_uF_per_km = C_F_per_m * 1e9

    # Step 5: Print the detailed calculation process and final result
    print("The formula for capacitance (C) per phase to screen is:")
    print("C = (2 * pi * epsilon_0 * epsilon_r) / ln( (R / r_w) * (1 - (m/R)^6)^(1/3) )")
    print("\n--- Calculation Breakdown ---")

    print(f"\nThe equation with the given values is:")
    print(f"C = (2 * {math.pi:.5f} * {epsilon_0:.5e} * {epsilon_r}) / ln( ({R} / {r_w}) * (1 - ({m}/{R})^6)^(1/3) )")
    
    # Print the evaluation of each part of the equation
    print(f"\nNumerator = 2 * pi * {epsilon:.5e} = {numerator:.5e} F/m")
    
    print(f"\nDenominator = ln( ({R_div_rw:.5f}) * (1 - ({m_div_R:.5f})^6)^(1/3) )")
    print(f"            = ln( {R_div_rw:.5f} * (1 - {m_div_R_pow6:.5f})^(1/3) )")
    print(f"            = ln( {R_div_rw:.5f} * ({term1_in_root:.5f})^(1/3) )")
    print(f"            = ln( {R_div_rw:.5f} * {term2_root:.5f} )")
    print(f"            = ln( {ln_argument:.5f} )")
    print(f"            = {denominator:.5f}")
            
    print("\n--- Final Result ---")
    print(f"C = {numerator:.5e} / {denominator:.5f} = {C_F_per_m:.5e} F/m")
    
    print(f"\nIn microfarads per kilometer, the capacitance is:")
    print(f"C = {C_F_per_m:.5e} F/m * 10^9 = {C_uF_per_km:.4f} µF/km")
    
    # Return the final numeric answer for parsing
    return C_uF_per_km

if __name__ == '__main__':
    final_capacitance = calculate_capacitance()
    # The final answer is wrapped for easy extraction.
    # print(f"\n<<<{final_capacitance:.4f}>>>")
    # To avoid double printing and confusion, the numeric-only answer is printed here.
    # The user can see the formatted answer above and use the raw value below if needed.

# Execute the function
calculate_capacitance()