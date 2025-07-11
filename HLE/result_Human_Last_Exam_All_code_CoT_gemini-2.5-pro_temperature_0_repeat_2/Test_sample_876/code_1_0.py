import math

def calculate_cable_capacitance():
    """
    Calculates the capacitance of a three-phase cable with a common screen.
    """
    # --- Given values ---
    # Radius of the wires in mm
    r_w_mm = 11.25
    # Distance from the center of a wire to the center of the cable in mm
    m_mm = 17.32
    # External radius of the cable (inner radius of the screen) in mm
    R_mm = 32.32
    # Relative permittivity of the insulator
    epsilon_r = 4.2

    # --- Constants ---
    # Permittivity of free space in F/m
    epsilon_0 = 8.854e-12

    # --- Unit Conversion ---
    # Convert all length units from mm to m for calculation
    r_w = r_w_mm / 1000
    m = m_mm / 1000
    R = R_mm / 1000

    # --- Calculation ---
    # The formula for capacitance per phase in F/m is:
    # C = (2 * pi * epsilon_0 * epsilon_r) / ln((R^3 - m^3) / (3 * r_w * m * R))

    # Calculate the argument of the natural logarithm
    log_arg_numerator = R**3 - m**3
    log_arg_denominator = 3 * r_w * m * R
    log_arg = log_arg_numerator / log_arg_denominator

    # Calculate the capacitance in Farads per meter (F/m)
    C_F_per_m = (2 * math.pi * epsilon_0 * epsilon_r) / math.log(log_arg)

    # Convert the result from F/m to microfarads per kilometer (uF/km)
    # Conversion factor: 1 F/m = 10^9 uF/km
    C_uF_per_km = C_F_per_m * 1e9

    # --- Output ---
    # Print the final equation with all the numbers substituted.
    # The calculation is for C in uF/km, so we include the 10^9 conversion factor.
    print("The capacitance C in uF/km is calculated using the formula:")
    print("C = (2 * pi * epsilon_0 * epsilon_r) / ln((R^3 - m^3) / (3 * r_w * m * R)) * 10^9")
    print("\nSubstituting the values:")
    print(f"C = (2 * {math.pi:.5f} * {epsilon_0} * {epsilon_r}) / ln((({R})^3 - ({m})^3) / (3 * {r_w} * {m} * {R})) * 10^9")
    print(f"\nThe calculated capacitance is:")
    print(f"C = {C_uF_per_km:.4f} uF/km")

if __name__ == '__main__':
    calculate_cable_capacitance()
<<<0.5660>>>