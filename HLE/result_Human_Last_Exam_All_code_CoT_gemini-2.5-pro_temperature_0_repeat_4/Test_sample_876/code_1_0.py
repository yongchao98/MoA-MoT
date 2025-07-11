import math

def calculate_cable_capacitance():
    """
    Calculates the capacitance of a three-phase cable with a common screen.
    """
    # 1. Define constants and given values
    # Permittivity of free space in Farads per meter (F/m)
    epsilon_0 = 8.854e-12
    
    # Given dimensions in millimeters (mm)
    r_w_mm = 11.25  # wire radius
    m_mm = 17.32    # distance from wire center to cable center
    R_mm = 32.32    # external radius of the cable (inner radius of the screen)
    
    # Relative permittivity of the insulator (dimensionless)
    epsilon_r = 4.2

    # 2. Convert dimensions from mm to meters (m) for SI unit consistency
    r_w = r_w_mm / 1000
    m = m_mm / 1000
    R = R_mm / 1000

    # 3. Calculate the absolute permittivity of the insulator
    epsilon = epsilon_r * epsilon_0

    # 4. Calculate capacitance in Farads per meter (F/m)
    # The formula approximates the capacitance per phase as that of an eccentric core.
    # C = (2 * pi * epsilon) / ln((R^2 - m^2) / (r_w * R))
    try:
        numerator = 2 * math.pi * epsilon
        log_argument = (R**2 - m**2) / (r_w * R)
        
        if log_argument <= 1:
            print("Error: Invalid geometry. The argument of the natural logarithm must be greater than 1.")
            return

        denominator = math.log(log_argument)
        C_F_per_m = numerator / denominator
    except (ValueError, ZeroDivisionError) as e:
        print(f"An error occurred during calculation: {e}")
        return

    # 5. Convert capacitance from F/m to microfarads per kilometer (µF/km)
    # Conversion factor: 1 F/m = 1e9 µF/km
    C_muF_per_km = C_F_per_m * 1e9

    # 6. Print the final equation with numbers and the result
    print("The final equation for capacitance (C) with the given values is:")
    # Using f-string to substitute numerical values into the equation string
    # The equation calculates the value in F/m, which is then multiplied by 10^9 for the final unit.
    print(f"C = [ (2 * {math.pi:.5f} * {epsilon_r} * {epsilon_0}) / ln(( {R}**2 - {m}**2) / ({r_w} * {R})) ] * 10^9 µF/km")
    
    print("\nResult:")
    print(f"The capacitance of the cable is {C_muF_per_km:.3f} µF/km.")

if __name__ == '__main__':
    calculate_cable_capacitance()