import math

def solve_triangle_area():
    """
    This function calculates and prints the formula for the area of triangle T(t).

    The derivation steps are as follows:
    1. The radius of the circumscribed circle R is 10, so the side length of the regular hexagon 's' is also 10.
    2. The rotation of the hexagon does not change the area of the inscribed triangle, so we can ignore it.
    3. The triangle T(t) starts as equilateral and remains equilateral due to the symmetry of the setup.
    4. We need to find the side length L(t) of the triangle T(t). We can set up a coordinate system for the hexagon.
       Let the vertices of the hexagon be at (10,0), (5, 5*sqrt(3)), (-5, 5*sqrt(3)), etc.
       The vertices of the triangle T(0) are midpoints of alternating sides, e.g., BC, DE, FA.
       - M1 on BC: (0, 5*sqrt(3))
       - M2 on DE: (-7.5, -2.5*sqrt(3))
       - M3 on FA: (7.5, -2.5*sqrt(3))
       The squared side length of T(0) is L(0)^2 = 15^2 = 225.
    5. The vertices of T(t) move from these midpoints with speed v=1. For example, the vertex P1(t) on side BC moves from M1 towards C. Its coordinates become P1(t) = (-t, 5*sqrt(3)).
       Similarly, we find the coordinates for the other two vertices.
    6. We calculate the squared distance between two vertices of T(t), for instance, P1(t) and P3(t).
       L(t)^2 = (x3(t) - x1(t))^2 + (y3(t) - y1(t))^2
       After algebraic simplification, this gives: L(t)^2 = 3*t^2 + 225.
    7. The area of an equilateral triangle is Area = (sqrt(3)/4) * L^2.
       Substituting L(t)^2, we get Area(t) = (sqrt(3)/4) * (3*t^2 + 225).
    """

    # Coefficients for the area formula Area(t) = a*t^2 + b
    # a = (3 * sqrt(3)) / 4
    # b = (225 * sqrt(3)) / 4
    a_coeff = (3 * math.sqrt(3)) / 4
    b_coeff = (225 * math.sqrt(3)) / 4

    # The final equation for the area of the triangle T(t)
    print("The area of the triangle T(t) as a function of time t is given by the formula:")
    print(f"Area(t) = (3 * sqrt(3) / 4) * t^2 + (225 * sqrt(3) / 4)")
    print("\nBreaking down the equation with calculated coefficients:")
    print(f"The coefficient for t^2 is: {a_coeff}")
    print(f"The constant term is: {b_coeff}")
    print(f"\nSo, the final equation is:")
    print(f"Area(t) = {a_coeff} * t^2 + {b_coeff}")
    print("\nThis can also be written in a factored form:")
    print("Area(t) = (sqrt(3)/4) * (3*t^2 + 225)")

solve_triangle_area()