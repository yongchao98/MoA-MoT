import math

def solve_playdough_gravity():
    """
    Calculates the maximum distance from point A to the surface of a 1 m^3
    volume shaped to create the strongest gravitational field at A.
    """
    
    # The derived formula for the maximum distance R_max is (15 / (4 * pi))^(1/3).
    # Here are the numbers from the formula:
    numerator = 15
    denominator_factor = 4
    pi = math.pi
    exponent_numerator = 1
    exponent_denominator = 3
    
    # Perform the calculation
    result = (numerator / (denominator_factor * pi)) ** (exponent_numerator / exponent_denominator)
    
    # Print the equation with all its numeric components, as requested.
    print("The formula for the furthest distance (R_max) is:")
    print(f"R_max = ({numerator} / ({denominator_factor} * {pi})) ^ ({exponent_numerator}/{exponent_denominator})")
    
    # Print the final calculated result.
    print("\nThe furthest distance is:")
    print(result)

solve_playdough_gravity()