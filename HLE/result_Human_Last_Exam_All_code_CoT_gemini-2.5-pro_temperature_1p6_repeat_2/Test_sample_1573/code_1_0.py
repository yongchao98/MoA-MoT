import math

def solve_chair_problem():
    """
    This script solves the chair problem by checking if the five leg endpoints
    are concyclic, which is a necessary condition for them to be cospherical.
    """
    # The five leg endpoint coordinates in the chair's reference plane.
    p1 = (0, 0)
    p2 = (2, 0)
    p3 = (2, 2)
    p4 = (0, 2)
    p5 = (1, 4)

    print("To solve this problem, we check if the five coplanar points of the leg tips are concyclic.")
    print("Let's find the circle defined by P1(0,0), P2(2,0), and P4(0,2).")
    print("The equation of a circle is (x - h)^2 + (y - k)^2 = r^2.\n")

    # From P1(0,0): h^2 + k^2 = r^2
    # From P2(2,0): (2-h)^2 + k^2 = r^2  => 4 - 4h + h^2 + k^2 = r^2 => 4 - 4h = 0 => h = 1
    # From P4(0,2): h^2 + (2-k)^2 = r^2  => h^2 + 4 - 4k + k^2 = r^2 => 4 - 4k = 0 => k = 1
    h = 1
    k = 1

    # Calculate r^2 using P1(0,0) and the center (h,k)
    r_squared = (p1[0] - h)**2 + (p1[1] - k)**2

    print(f"The center of the circle is (h,k) = ({h},{k}).")
    print(f"The radius squared is r^2 = ({p1[0]} - {h})^2 + ({p1[1]} - {k})^2 = {r_squared}.")
    print(f"So, the equation of the circle is (x - {h})^2 + (y - {k})^2 = {int(r_squared)}.\n")

    print("Now, we check if the other points, P3 and P5, are on this circle.")

    # Check P3(2,2)
    dist_sq_p3 = (p3[0] - h)**2 + (p3[1] - k)**2
    print(f"\nChecking point P3{p3}:")
    print(f"  ({p3[0]} - {h})^2 + ({p3[1]} - {k})^2 = {int(dist_sq_p3)}")
    if math.isclose(dist_sq_p3, r_squared):
        print(f"  Result: {int(dist_sq_p3)} == {int(r_squared)}. Point P3 is on the circle.")
    else:
        print(f"  Result: {int(dist_sq_p3)} != {int(r_squared)}. Point P3 is NOT on the circle.")

    # Check P5(1,4)
    dist_sq_p5 = (p5[0] - h)**2 + (p5[1] - k)**2
    print(f"\nChecking point P5{p5}:")
    print(f"  ({p5[0]} - {h})^2 + ({p5[1]} - {k})^2 = {int(dist_sq_p5)}")
    if math.isclose(dist_sq_p5, r_squared):
        print(f"  Result: {int(dist_sq_p5)} == {int(r_squared)}. Point P5 is on the circle.")
    else:
        print(f"  Result: {int(dist_sq_p5)} != {int(r_squared)}. Point P5 is NOT on the circle.")

    print("\nConclusion: Since Point P5 does not lie on the same circle as P1, P2, P3, and P4, the five points are not concyclic.")
    print("Therefore, they cannot be cospherical. It is impossible for all five legs to touch a spherical surface simultaneously.")
    print("The number of such locations is 0.")

solve_chair_problem()