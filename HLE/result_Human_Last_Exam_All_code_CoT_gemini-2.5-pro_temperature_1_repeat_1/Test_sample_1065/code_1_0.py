import math

def calculate_center_of_mass():
    """
    This function calculates the normalized horizontal and vertical coordinates 
    of the center of mass of the candy string.
    
    The formulas derived are:
    X_cm / R = (6 + 3*sqrt(3)) / (4*pi + 6)
    Z_cm / R = (12*sqrt(3) + 3) / (8*pi + 12)
    """

    # --- Calculation for the horizontal coordinate (X_cm/R) ---
    
    # Numerator parts for X: 6 + 3*sqrt(3)
    num_x_part1 = 6.0
    num_x_part2 = 3 * math.sqrt(3)
    numerator_x = num_x_part1 + num_x_part2

    # Denominator parts for X: 4*pi + 6
    den_x_part1 = 4 * math.pi
    den_x_part2 = 6.0
    denominator_x = den_x_part1 + den_x_part2
    
    x_cm_norm = numerator_x / denominator_x

    print("--- Calculating the horizontal coordinate (X_cm/R) ---")
    print("Equation for X_cm/R: (A + B) / (C + D)")
    print(f"The number for A is: {num_x_part1}")
    print(f"The number for B (3*sqrt(3)) is: {num_x_part2}")
    print(f"The number for C (4*pi) is: {den_x_part1}")
    print(f"The number for D is: {den_x_part2}")
    print(f"Resulting X_cm/R = ({numerator_x}) / ({denominator_x}) = {x_cm_norm}\n")

    # --- Calculation for the vertical coordinate (Z_cm/R) ---

    # Numerator parts for Z: 12*sqrt(3) + 3
    num_z_part1 = 12 * math.sqrt(3)
    num_z_part2 = 3.0
    numerator_z = num_z_part1 + num_z_part2

    # Denominator parts for Z: 8*pi + 12
    den_z_part1 = 8 * math.pi
    den_z_part2 = 12.0
    denominator_z = den_z_part1 + den_z_part2

    z_cm_norm = numerator_z / denominator_z

    print("--- Calculating the vertical coordinate (Z_cm/R) ---")
    print("Equation for Z_cm/R: (E + F) / (G + H)")
    print(f"The number for E (12*sqrt(3)) is: {num_z_part1}")
    print(f"The number for F is: {num_z_part2}")
    print(f"The number for G (8*pi) is: {den_z_part1}")
    print(f"The number for H is: {den_z_part2}")
    print(f"Resulting Z_cm/R = ({numerator_z}) / ({denominator_z}) = {z_cm_norm}\n")

    # --- Final Answer ---
    print("The raw numbers for the horizontal and vertical coordinates, separated by a comma:")
    print(f"{x_cm_norm},{z_cm_norm}")

if __name__ == '__main__':
    calculate_center_of_mass()