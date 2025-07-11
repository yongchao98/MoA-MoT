import math

def calculate_cable_capacitance():
    """
    Calculates the capacitance of a three-phase cable with a common screen
    and prints the detailed steps and final result.
    """
    # --- Given values ---
    # Wire radius in mm
    r_w_mm = 11.25
    # Distance from wire center to cable center in mm
    m_mm = 17.32
    # External radius of the cable in mm
    R_mm = 32.32
    # Relative permittivity of the insulator
    epsilon_r = 4.2

    # --- Constants ---
    # Permittivity of free space in F/m
    epsilon_0 = 8.854e-12
    # Constant pi
    pi = math.pi

    # --- Convert units to SI (meters) ---
    r_w = r_w_mm / 1000.0
    m = m_mm / 1000.0
    R = R_mm / 1000.0

    # --- Calculation Steps ---
    # Formula: C = (2 * pi * epsilon_0 * epsilon_r) / ln((R^2 - m^2) / (r_w * R))

    # Step 1: Calculate the numerator
    numerator = 2 * pi * epsilon_0 * epsilon_r

    # Step 2: Calculate the argument of the natural logarithm
    # Using mm for this ratio is fine as units cancel, but we'll use SI units for consistency.
    log_argument = (R**2 - m**2) / (r_w * R)

    # Step 3: Calculate the denominator
    denominator = math.log(log_argument)

    # Step 4: Calculate capacitance in Farads per meter (F/m)
    C_F_per_m = numerator / denominator

    # Step 5: Convert capacitance to microfarads per kilometer (uF/km)
    # Conversion factor: 1 F/m = 1e9 uF/km
    C_muF_per_km = C_F_per_m * 1e9

    # --- Output Results ---
    print("This script calculates the capacitance of a three-phase cable with a common screen.")
    print("\nThe formula used is: C = (2 * pi * \u03B5\u2080 * \u03B5\u1D63) / ln((R\u00B2 - m\u00B2) / (r\u209b * R))")
    
    print("\n--- Equation with numerical values ---")
    print(f"C = (2 * {pi:.5f} * {epsilon_0:.4e} F/m * {epsilon_r}) / ln((({R} m)\u00B2 - ({m} m)\u00B2) / (({r_w} m) * ({R} m)))")

    print("\n--- Calculation Breakdown ---")
    print(f"Numerator = 2 * {pi:.5f} * {epsilon_0:.4e} * {epsilon_r} = {numerator:.4e} F/m")
    # To show the ratio calculation clearly, we can show the calculation with the mm values
    print(f"Logarithm Argument = ({R_mm}\u00B2 - {m_mm}\u00B2) / ({r_w_mm} * {R_mm}) = {log_argument:.4f}")
    print(f"Denominator = ln({log_argument:.4f}) = {denominator:.4f}")
    print(f"Capacitance (C) = {numerator:.4e} / {denominator:.4f} = {C_F_per_m:.4e} F/m")

    print("\n--- Final Result ---")
    print(f"The capacitance of the cable is {C_muF_per_km:.3f} \u03BCF/km.")

if __name__ == '__main__':
    calculate_cable_capacitance()
<<<0.326>>>