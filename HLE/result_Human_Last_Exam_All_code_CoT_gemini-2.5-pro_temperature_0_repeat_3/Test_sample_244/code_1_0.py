def solve_knot_proportion():
    """
    This function calculates the proportion of 7-crossing knots that are hyperbolic.
    """
    # Step 1: Determine the total number of 7-crossing knots.
    # From knot theory, there are 7 distinct prime knots with a crossing number of 7.
    # They are denoted 7_1, 7_2, 7_3, 7_4, 7_5, 7_6, and 7_7.
    total_knots = 7

    # Step 2 & 3: Identify non-hyperbolic 7-crossing knots.
    # A prime knot is non-hyperbolic if it is a torus knot.
    # We check for torus knots T(p,q) with crossing number 7.
    # The crossing number of T(p,q) is min(p(q-1), q(p-1)).
    # For the crossing number to be 7 (a prime number), the only solution for integers p,q > 1
    # is the T(2,7) knot, which has crossing number min(2(7-1), 7(2-1)) = min(12, 7) = 7.
    # This knot is conventionally known as 7_1.
    # There are no satellite knots with 7 crossings.
    non_hyperbolic_knots = 1 # This is the 7_1 knot.

    # Step 4: Calculate the number of hyperbolic knots.
    hyperbolic_knots = total_knots - non_hyperbolic_knots

    # Step 5: Calculate the proportion.
    proportion = hyperbolic_knots / total_knots

    # Step 6: Print the results, including the final equation.
    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of non-hyperbolic knots (torus knots): {non_hyperbolic_knots}")
    print(f"Number of hyperbolic knots: {hyperbolic_knots}")
    print("\nThe proportion of hyperbolic knots is calculated as:")
    print(f"Proportion = (Number of hyperbolic knots) / (Total number of knots)")
    print(f"Proportion = {hyperbolic_knots} / {total_knots}")
    print(f"The decimal value is approximately: {proportion:.4f}")

solve_knot_proportion()