import math

def calculate_cable_capacitance():
    """
    Calculates the capacitance of a three-phase cable with a common screen and
    prints the steps and the final equation.
    """
    # Step 1: Define given values and constants
    r_w_mm = 11.25  # Wire radius in mm
    m_mm = 17.32    # Distance from wire center to cable center in mm
    R_mm = 32.32    # External radius of the cable in mm
    epsilon_r = 4.2 # Relative permittivity of the insulator
    epsilon_0 = 8.854e-12 # Permittivity of free space in F/m

    # Step 2: Convert all lengths to SI units (meters) for calculation
    r_w = r_w_mm / 1000
    m = m_mm / 1000
    R = R_mm / 1000

    # Step 3: Calculate the capacitance using the formula
    # C [F/m] = (2 * pi * epsilon_0 * epsilon_r) / ln((R^2 - m^2) / (r_w * R))
    
    # Calculate terms for clarity
    numerator = 2 * math.pi * epsilon_0 * epsilon_r
    log_arg_num = R**2 - m**2
    log_arg_den = r_w * R
    log_argument = log_arg_num / log_arg_den
    denominator = math.log(log_argument)
    
    # Capacitance in Farads per meter (F/m)
    C_F_per_m = numerator / denominator

    # Step 4: Convert the result to microfarads per kilometer (μF/km)
    # Conversion factor: 1 F/m = 10^9 μF/km
    C_uF_per_km = C_F_per_m * 1e9

    # Step 5: Print the detailed calculation and result
    print("--- Final Equation and Calculation ---")
    print("The formula for capacitance per phase is: C = (2 * π * ε₀ * εᵣ) / ln((R² - m²) / (r_w * R))\n")
    
    print("Plugging in the values (with lengths in meters):")
    # Show the formula with all numbers substituted
    print(f"C (μF/km) = [ (2 * {math.pi:.5f} * {epsilon_0:.3e} * {epsilon_r}) / ln( (({R})² - ({m})²) / ({r_w} * {R}) ) ] * 10^9")
    
    # Show the evaluation step-by-step
    print(f"C (μF/km) = [ ({numerator:.5e}) / ln( ({log_arg_num:.5e}) / ({log_arg_den:.5e}) ) ] * 10^9")
    print(f"C (μF/km) = [ ({numerator:.5e}) / ln({log_argument:.5f}) ] * 10^9")
    print(f"C (μF/km) = [ ({numerator:.5e}) / {denominator:.5f} ] * 10^9")
    print(f"C (μF/km) = [ {C_F_per_m:.5e} ] * 10^9\n")

    print(f"The final capacitance of the cable is {C_uF_per_km:.3f} μF/km.")
    
    return C_uF_per_km

# Execute the function to get the answer
final_answer = calculate_cable_capacitance()
print(f"<<<{final_answer:.3f}>>>")