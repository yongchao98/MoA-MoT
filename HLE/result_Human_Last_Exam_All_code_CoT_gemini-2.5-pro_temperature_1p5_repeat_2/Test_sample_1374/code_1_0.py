import math

def solve_playdough_problem():
    """
    Calculates the maximum distance from point A for a volume of playdough
    shaped to create the strongest gravitational field at A.

    The shape that maximizes the gravitational field at its apex (point A)
    is a solid of revolution whose surface is described by r = k * sqrt(cos(θ)).
    The volume of this shape is V = (4 * pi * k^3) / 15.

    Given V = 1 cubic meter, we solve for k, which represents the maximum
    distance from A to the surface of the playdough.
    """
    
    # Constants from the volume formula V = (4 * pi * k^3) / 15
    # We solve for k when V = 1: k = (15 / (4 * pi))^(1/3)
    numerator = 15
    denominator_factor = 4
    exponent_numerator = 1
    exponent_denominator = 3
    
    print("The shape that maximizes the gravitational field at a point A is a solid of revolution.")
    print("Its surface is described by r = k * sqrt(cos(θ)), and its volume is V = (4 * pi * k^3) / 15.")
    print("For a volume of 1 cubic meter, we can solve for k.")
    print("The furthest point on the surface is at θ=0, where r_max = k.")
    print("\nThe equation for this distance is:")
    print(f"r_max = k = ({numerator} / ({denominator_factor} * pi)) ^ ({exponent_numerator}/{exponent_denominator})")

    # Perform the calculation
    k = (numerator / (denominator_factor * math.pi)) ** (exponent_numerator / exponent_denominator)

    print("\nCalculating the result:")
    print(f"The furthest point on the surface of the playdough is {k:.6f} meters from point A.")

solve_playdough_problem()