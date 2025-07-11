import math

def solve_gravity_shape_problem():
    """
    Calculates the maximum distance from a point A to the surface of a 1 m^3 volume
    shaped to maximize the gravitational field at A.
    """
    
    # The shape is defined by r = a * sqrt(cos(θ)). Its volume V is (4*pi*a^3)/15.
    # Given V = 1, we solve for 'a', which is the maximum distance.
    # a = (15 * V / (4 * pi))^(1/3)
    
    # Define the numbers in the final equation
    volume = 1.0
    numerator = 15.0
    denominator_factor1 = 4.0
    denominator_factor2 = math.pi
    exponent_numerator = 1.0
    exponent_denominator = 3.0
    
    # Calculate the result
    max_distance = ( (numerator * volume) / (denominator_factor1 * denominator_factor2) ) ** (exponent_numerator / exponent_denominator)

    # Print the final equation with its components
    print(f"The equation for the furthest distance 'd' is derived from the shape's volume.")
    print(f"d = ({numerator} / ({denominator_factor1} * {denominator_factor2:.6f})) ^ ({exponent_numerator:.0f}/{exponent_denominator:.0f})")
    
    # Print the final result
    print(f"\nThe furthest point on the surface of the playdough is {max_distance:.6f} meters from point A.")

solve_gravity_shape_problem()