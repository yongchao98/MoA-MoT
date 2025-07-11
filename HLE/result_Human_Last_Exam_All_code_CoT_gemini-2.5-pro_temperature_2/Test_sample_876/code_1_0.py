import math

def calculate_cable_capacitance():
    """
    Calculates the capacitance of a three-phase cable with a common screen.
    """
    # --- Given values in mm ---
    r_w_mm = 11.25  # radius of the wire in mm
    m_mm = 17.32    # distance from wire center to cable center in mm
    R_mm = 32.32    # external radius of the cable (inner screen radius) in mm
    epsilon_r = 4.2 # relative permittivity of the insulator

    # --- Constants ---
    epsilon_0 = 8.854e-12  # Permittivity of free space in F/m

    # --- Convert units to SI (meters) for calculation ---
    r_w = r_w_mm / 1000
    m = m_mm / 1000
    R = R_mm / 1000

    # --- Step-by-step calculation ---

    # 1. Calculate the permittivity of the insulator
    epsilon = epsilon_r * epsilon_0

    # 2. Calculate the term inside the logarithm
    log_argument_numerator = R**3 - m**3
    log_argument_denominator = 3 * r_w * m**2
    log_argument = log_argument_numerator / log_argument_denominator

    # 3. Calculate the natural logarithm of the argument
    ln_val = math.log(log_argument)

    # 4. Calculate the capacitance in Farads per meter (F/m)
    C_F_per_m = (2 * math.pi * epsilon) / ln_val

    # 5. Convert the capacitance to microfarads per kilometer (uF/km)
    C_uF_per_km = C_F_per_m * 1e9

    # --- Print the results ---
    print("This script calculates the per-phase capacitance of a three-phase cable with a common screen.")
    print("The formula is: C = (2 * pi * epsilon) / ln((R^3 - m^3) / (3 * r_w * m^2))\n")
    print("--- Calculation Steps ---")
    
    # Print each number used in the calculation as requested
    print("1. Permittivity (epsilon):")
    print(f"   epsilon = epsilon_r * epsilon_0 = {epsilon_r} * {epsilon_0:.4e} F/m = {epsilon:.4e} F/m")

    print("\n2. Argument of the natural logarithm (ln):")
    # Using mm values for clearer representation of the inputs
    log_arg_num_mm = R_mm**3 - m_mm**3
    log_arg_den_mm = 3 * r_w_mm * m_mm**2
    print(f"   Argument = (R^3 - m^3) / (3 * r_w * m^2)")
    print(f"            = ({R_mm:.2f}^3 - {m_mm:.2f}^3) / (3 * {r_w_mm:.2f} * {m_mm:.2f}^2)")
    print(f"            = ({R_mm**3:.2f} - {m_mm**3:.2f}) / ({log_arg_den_mm:.2f})")
    print(f"            = {log_argument:.4f}")

    print("\n3. Capacitance in F/m:")
    print(f"   C = (2 * pi * {epsilon:.4e}) / ln({log_argument:.4f})")
    print(f"     = {2 * math.pi * epsilon:.4e} / {ln_val:.4f} = {C_F_per_m:.4e} F/m")

    print("\n--- Final Result ---")
    print("To get the result in uF/km, we multiply by 10^9:")
    print(f"Capacitance = {C_F_per_m:.4e} F/m * 1e9 = {C_uF_per_km:.4f} uF/km")

    # Return the final numeric value for the answer tag
    return C_uF_per_km

# Execute the calculation and print the final answer tag
final_capacitance = calculate_cable_capacitance()
print(f"\n<<<C = {final_capacitance:.4f} uF/km>>>")