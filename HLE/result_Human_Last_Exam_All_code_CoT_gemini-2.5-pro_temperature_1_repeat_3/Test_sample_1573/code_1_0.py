import math

def solve():
    """
    This script determines if the five leg endpoints of a chair can all touch a spherical surface.
    This is possible if and only if the five coplanar points are concyclic (lie on the same circle).
    """

    # The coordinates of the five leg endpoints in a plane
    p1 = (0, 0)
    p2 = (2, 0)
    p3 = (2, 2)
    p4 = (0, 2)
    p5 = (1, 4)

    print("Step 1: Define the five leg endpoint coordinates.")
    print(f"P1 = {p1}, P2 = {p2}, P3 = {p3}, P4 = {p4}, P5 = {p5}\n")

    print("Step 2: Determine the equation of a circle from three points (P1, P2, P4).")
    # The center of the circle (h, k) is the intersection of the perpendicular bisectors.
    # Perpendicular bisector of P1(0,0) and P2(2,0) is x = 1.
    # Perpendicular bisector of P1(0,0) and P4(0,2) is y = 1.
    h, k = 1, 1
    print(f"The center of the circle is (h, k) = ({h}, {k}).")

    # The radius squared (r^2) is the squared distance from the center to any of the points.
    # Using P1(0,0):
    r_squared = (p1[0] - h)**2 + (p1[1] - k)**2
    print(f"The radius squared is r^2 = ({p1[0]} - {h})^2 + ({p1[1]} - {k})^2 = {r_squared}.\n")
    
    print(f"The equation of the circle is (x - {h})^2 + (y - {k})^2 = {r_squared}.\n")

    print("Step 3: Check if the other points (P3 and P5) lie on this circle.")

    # Check P3(2,2)
    dist_sq_p3 = (p3[0] - h)**2 + (p3[1] - k)**2
    print(f"For P3{p3}: ({p3[0]} - {h})^2 + ({p3[1]} - {k})^2 = {dist_sq_p3}")
    if math.isclose(dist_sq_p3, r_squared):
        print("Result: P3 lies on the circle.\n")
    else:
        print("Result: P3 does NOT lie on the circle.\n")

    # Check P5(1,4)
    dist_sq_p5 = (p5[0] - h)**2 + (p5[1] - k)**2
    print(f"For P5{p5}: ({p5[0]} - {h})^2 + ({p5[1]} - {k})^2 = {dist_sq_p5}")
    if math.isclose(dist_sq_p5, r_squared):
        print("Result: P5 lies on the circle.\n")
    else:
        print(f"Result: P5 does NOT lie on the circle, since {dist_sq_p5} != {r_squared}.\n")

    print("Conclusion:")
    print("The five points are coplanar but not concyclic.")
    print("Therefore, they cannot all lie on the surface of a sphere simultaneously.")
    print("The number of locations where all five legs can touch the surface is 0.")

solve()