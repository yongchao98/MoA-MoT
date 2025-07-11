import math

def solve_pumpkin_string_cm():
    """
    Calculates the horizontal and vertical coordinates of the center of mass
    of the candy string as described in the problem.
    """
    
    # We are asked for the raw numbers, so we can assume R=1.
    R = 1.0
    
    # Define mathematical constants
    pi = math.pi
    sqrt3 = math.sqrt(3)
    
    # The final equations for the coordinates of the center of mass are derived as:
    # r_cm = (3 * R * (2 + sqrt(3))) / (2 * (2 * pi + 3))
    # z_cm = (3 * R * (4 * sqrt(3) + 1)) / (4 * (2 * pi + 3))
    
    # Let's show the equation and its components for clarity
    print("The horizontal coordinate (r_cm) is calculated by the formula:")
    print("r_cm = (3 * (2 + sqrt(3))) / (2 * (2 * pi + 3))")
    
    # Numerator and denominator for r_cm
    r_num = 3 * (2 + sqrt3)
    r_den = 2 * (2 * pi + 3)
    r_cm = r_num / r_den
    
    print(f"Numerator value: {r_num}")
    print(f"Denominator value: {r_den}")
    
    print("\nThe vertical coordinate (z_cm) is calculated by the formula:")
    print("z_cm = (3 * (4 * sqrt(3) + 1)) / (4 * (2 * pi + 3))")
    
    # Numerator and denominator for z_cm
    z_num = 3 * (4 * sqrt3 + 1)
    z_den = 4 * (2 * pi + 3)
    z_cm = z_num / z_den
    
    print(f"Numerator value: {z_num}")
    print(f"Denominator value: {z_den}")
    
    # Print the final raw numbers separated by a comma
    print("\nThe raw horizontal and vertical coordinates are:")
    print(f"{r_cm},{z_cm}")

solve_pumpkin_string_cm()