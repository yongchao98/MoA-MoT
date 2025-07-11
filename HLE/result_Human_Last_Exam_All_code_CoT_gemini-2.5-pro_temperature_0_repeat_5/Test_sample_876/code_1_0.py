import math

def calculate_capacitance():
    """
    Calculates the capacitance of a three-phase cable with a common screen.
    """
    # Step 1: Define parameters and constants
    # Given values
    r_w_mm = 11.25  # wire radius in mm
    m_mm = 17.32   # distance from wire center to cable center in mm
    R_mm = 32.32   # inner radius of the common screen in mm
    epsilon_r = 4.2 # relative permittivity of the insulator

    # Constants
    epsilon_0 = 8.854e-12 # permittivity of free space in F/m

    # Convert units from mm to m for calculation
    r_w = r_w_mm / 1000
    m = m_mm / 1000
    R = R_mm / 1000

    # Step 2: Calculate terms for the formula
    # The formula is C = (2 * pi * epsilon_0 * epsilon_r) / ln((R^2 - m^2) / (r_w * R))
    epsilon = epsilon_r * epsilon_0
    numerator_C = 2 * math.pi * epsilon

    log_arg_numerator = R**2 - m**2
    log_arg_denominator = r_w * R
    log_arg = log_arg_numerator / log_arg_denominator
    denominator_C = math.log(log_arg)

    # Step 3: Calculate capacitance in F/m
    C_F_per_m = numerator_C / denominator_C

    # Step 4: Convert to μF/km
    # 1 F/m = 10^9 μF/km
    C_muF_per_km = C_F_per_m * 1e9

    # Step 5: Print the calculation and result
    print("Calculation of the capacitance of a three-phase cable with a common screen.")
    print("-" * 70)
    print("Formula: C = (2 * pi * epsilon_0 * epsilon_r) / ln((R^2 - m^2) / (r_w * R))")
    print("\nGiven values:")
    print(f"Wire radius (r_w) = {r_w_mm} mm")
    print(f"Distance to center (m) = {m_mm} mm")
    print(f"Screen inner radius (R) = {R_mm} mm")
    print(f"Relative permittivity (epsilon_r) = {epsilon_r}")
    print(f"Permittivity of free space (epsilon_0) = {epsilon_0:.4e} F/m")
    print("-" * 70)

    print("Substituting the values into the formula (in SI units):")
    print(f"C = (2 * {math.pi:.4f} * {epsilon_0:.4e} * {epsilon_r}) / ln(({R:.5f}^2 - {m:.5f}^2) / ({r_w:.5f} * {R:.5f}))")
    print(f"C = ({numerator_C:.4e}) / ln(({R**2:.6f} - {m**2:.6f}) / {log_arg_denominator:.6f})")
    print(f"C = ({numerator_C:.4e}) / ln({log_arg_numerator:.6f} / {log_arg_denominator:.6f})")
    print(f"C = ({numerator_C:.4e}) / ln({log_arg:.4f})")
    print(f"C = ({numerator_C:.4e}) / {denominator_C:.4f}")
    print(f"C = {C_F_per_m:.4e} F/m")
    
    print("\nConverting to microfarads per kilometer (μF/km):")
    print(f"C_μF/km = {C_F_per_m:.4e} F/m * 10^9 (μF/km)/(F/m)")
    print(f"C_μF/km = {C_muF_per_km:.3f} μF/km")
    print("-" * 70)
    print(f"The final capacitance is {C_muF_per_km:.3f} μF/km.")

if __name__ == '__main__':
    calculate_capacitance()
<<<0.326>>>