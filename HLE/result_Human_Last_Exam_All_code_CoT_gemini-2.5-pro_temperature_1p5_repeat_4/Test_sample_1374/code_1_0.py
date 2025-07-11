import math

def solve_playdough_problem():
    """
    Calculates the maximum distance from point A to the surface of a playdough shape
    that maximizes the gravitational field at A.

    The problem is solved by determining the optimal shape and then calculating its dimensions.
    The optimal shape that maximizes the gravitational field at a point A (the origin)
    is bounded by the surface r^2 = k * cos(theta).

    The volume (V) of this shape is V = (4 * pi * k^(3/2)) / 15.
    Given V = 1 m^3, we can solve for k.

    The furthest point on the surface occurs at theta = 0, where r_max = sqrt(k).
    This leads to the final equation for the maximum distance:
    r_max = (15 / (4 * pi))^(1/3)
    """

    # Define the constants from the final equation
    numerator = 15
    denominator_coeff = 4
    pi_val = math.pi
    exponent = 1/3

    # The final equation for the maximum distance r_max
    print("The final equation for the maximum distance r_max is:")
    print(f"r_max = ({numerator} / ({denominator_coeff} * pi))^(1/3)")
    print("\nPlugging in the values:")
    print(f"r_max = ({float(numerator)} / ({float(denominator_coeff)} * {pi_val}))^({exponent})")

    # Calculate the result
    r_max = (numerator / (denominator_coeff * pi_val))**exponent

    print("\nThe furthest distance from A (in meters) is:")
    print(r_max)

solve_playdough_problem()