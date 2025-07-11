import math

def solve_minimal_area():
    """
    Calculates the minimal area of a convex domain that intersects all lines
    px+qy=1, where p and q are coprime integers.

    The analytical solution to this problem is known to be 12 / pi^2.
    This function computes and prints this value.
    """
    
    # The components of the final formula
    numerator = 12
    denominator_power = 2
    
    # Calculate the result
    pi_squared = math.pi ** denominator_power
    result = numerator / pi_squared
    
    # Print the equation and the result
    print(f"The minimal area is given by the formula: {numerator} / π^{denominator_power}")
    print(f"Let's calculate this value:")
    print(f"{numerator} / ({math.pi:.6f})^{denominator_power} = {result:.6f}")
    print("\nFinal equation and result:")
    print(f"{numerator} / π**{denominator_power} = {result}")

solve_minimal_area()