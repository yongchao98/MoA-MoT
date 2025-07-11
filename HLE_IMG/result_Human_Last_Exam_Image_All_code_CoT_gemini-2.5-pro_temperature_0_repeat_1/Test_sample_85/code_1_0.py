import math

def solve_cone_distance():
    """
    Calculates the furthest distance from a point P on a cone's base circumference.

    The cone's properties are:
    - Base diameter = d
    - Slant height (distance from apex to P) = d

    The solution involves unrolling the cone and finding the maximum distance
    on the resulting 2D shape.
    """

    # The problem is symbolic, with the answer in terms of 'd'.
    # We will calculate the numerical coefficient for the final expression.

    # As derived in the plan, the unrolled cone is a semicircle of radius d.
    # The starting point P is on the arc, and the furthest point Q is at a corner.
    # In a coordinate system with the apex at (0,0), P is at (0,d) and Q is at (d,0).
    # The distance is the hypotenuse of a right triangle with both legs of length d.
    # Distance^2 = d^2 + d^2 = 2 * d^2
    # Distance = sqrt(2 * d^2) = d * sqrt(2)

    # The number in the equation is 2, inside the square root.
    number_in_equation = 2
    coefficient = math.sqrt(number_in_equation)

    print("The furthest distance on the cone's surface is found by calculating the longest straight-line distance on its unrolled 2D shape (a semicircle of radius d).")
    print("This distance forms the hypotenuse of a right triangle with two sides of length d.")
    print("\nThe calculation is as follows:")
    print("Distance^2 = d^2 + d^2")
    print(f"Distance^2 = {number_in_equation} * d^2")
    print(f"Distance = d * sqrt({number_in_equation})")
    print("\nTherefore, the furthest distance from point P is:")
    print(f"d * {coefficient}")

solve_cone_distance()