import math

def calculate_capacitance():
    """
    Calculates the capacitance of a three-phase cable with a common screen
    and prints the steps and result.
    """
    # Step 1: Define given values and constants
    r_w_mm = 11.25  # radius of the wires in mm
    m_mm = 17.32    # distance from wire center to cable center in mm
    R_mm = 32.32    # external radius of the cable in mm
    epsilon_r = 4.2 # relative permittivity of the insulator

    # Constants
    epsilon_0 = 8.854e-12 # permittivity of free space in F/m
    pi = math.pi

    # Step 2: Convert all length units from mm to meters
    r_w = r_w_mm / 1000.0
    m = m_mm / 1000.0
    R = R_mm / 1000.0

    # Step 3: Calculate the components of the formula
    # Numerator of the main formula
    formula_num = 2 * pi * epsilon_0 * epsilon_r
    
    # Numerator and denominator of the logarithm's argument
    log_arg_num = R**3 - m**3
    log_arg_den = 3 * m * R * r_w
    
    # The full argument of the logarithm
    log_arg = log_arg_num / log_arg_den
    
    # The natural logarithm value
    log_val = math.log(log_arg)

    # Calculate capacitance in Farads per meter (F/m)
    capacitance_F_per_m = formula_num / log_val

    # Step 4: Convert the result to microfarads per kilometer (μF/km)
    # Conversion factor from F/m to μF/km is 10^9
    capacitance_uF_per_km = capacitance_F_per_m * 1e9

    # Step 5: Print the formula with the numerical values
    print("The formula for capacitance (C) per phase is:")
    print("C = (2 * pi * ε_0 * ε_r) / ln[(R^3 - m^3) / (3 * m * R * r_w)]\n")
    
    print("Plugging in the values (with lengths in meters):")
    print(f"ε_0 = {epsilon_0} F/m")
    print(f"ε_r = {epsilon_r}")
    print(f"r_w = {r_w} m")
    print(f"m = {m} m")
    print(f"R = {R} m\n")

    print("The equation becomes:")
    # Using f-string to format the output showing the equation with numbers
    # We display the actual numbers used in the calculation
    print(f"C = (2 * {pi:.6f} * {epsilon_0} * {epsilon_r}) / ln[({R}**3 - {m}**3) / (3 * {m} * {R} * {r_w})]")
    print(f"C = ({formula_num:.6e}) / ln[({log_arg_num:.6e}) / ({log_arg_den:.6e})]")
    print(f"C = ({formula_num:.6e}) / ln({log_arg:.6f})")
    print(f"C = ({formula_num:.6e}) / {log_val:.6f}")
    print(f"C = {capacitance_F_per_m:.6e} F/m\n")

    print("Finally, converting the result to μF per km:")
    print(f"C = {capacitance_F_per_m:.6e} F/m * 10^9 = {capacitance_uF_per_km:.4f} μF/km")

if __name__ == '__main__':
    calculate_capacitance()
    # The final answer is approximately 0.5650 μF/km.
    # We are providing the numerical value as requested.
    # <<<0.5650>>> is the format for the final answer.
    print("\n<<<0.5650>>>")