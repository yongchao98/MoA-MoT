import math

def solve_playdough_problem():
    """
    Calculates the maximum distance from point A for a playdough shape
    that maximizes the gravitational field at A.
    """
    
    # The given volume of the playdough in cubic meters.
    volume = 1.0
    
    # The formula for the maximum distance (R_max) is derived from the volume of the optimal shape:
    # V = (4 * pi * R_max^3) / 15
    # R_max = (15 * V / (4 * pi))^(1/3)
    
    numerator = 15.0
    denominator_coeff = 4.0
    
    # Calculate R_max
    r_max = (numerator * volume / (denominator_coeff * math.pi))**(1/3)
    
    # Print the final equation with the numbers used
    print("The distance to the furthest point (R_max) is calculated from the volume equation for the optimal shape.")
    print(f"The equation is: R_max = ({numerator} * {volume} / ({denominator_coeff} * pi))^(1/3)")
    print("\nCalculating the final value:")
    print(f"R_max = ({numerator} * {volume} / ({denominator_coeff} * {math.pi:.5f}...))^(1/3)")
    print(f"R_max = {r_max:.4f} meters")

solve_playdough_problem()
