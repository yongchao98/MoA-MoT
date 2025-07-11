import math

def solve_radius():
    """
    Calculates the radius of a circle packing 14 smaller circles of radius 1,
    based on a (1, 5, 8) configuration.

    The formula for the radius R is derived as R = 3 + csc(pi/5).
    """
    
    # Numbers in the equation R = 3 + 1 / sin(pi/5)
    num1 = 3
    num2 = 1
    num3 = 5
    
    print(f"The calculation is based on the formula: R = {num1} + {num2} / sin(pi/{num3})")
    
    # Calculate csc(pi/5) which is 1 / sin(pi/5)
    csc_pi_over_5 = 1 / math.sin(math.pi / num3)
    
    # Calculate the final radius
    radius = num1 + csc_pi_over_5
    
    # Output the radius rounded to 4 significant digits
    print(f"The radius of the circle is: {radius:.4g}")

solve_radius()