def solve_knot_proportion():
    """
    This function calculates the proportion of hyperbolic knots among those
    with a crossing number of 7.
    """
    # 1. Total number of knot types with 7 crossings.
    # The knots are 7_1, 7_2, 7_3, 7_4, 7_5, 7_6, 7_7.
    total_knots = 7

    # 2. Identify non-hyperbolic knots.
    # All 7-crossing knots are prime. A prime knot is either a torus knot or a hyperbolic knot.
    # We identify the torus knot(s). The knot 7_1 is the T(2,7) torus knot.
    # It is the only torus knot with 7 crossings.
    num_non_hyperbolic = 1

    # 3. Calculate the number of hyperbolic knots.
    num_hyperbolic = total_knots - num_non_hyperbolic

    # 4. Calculate the proportion.
    proportion = num_hyperbolic / total_knots

    # Print the breakdown and the final equation.
    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of non-hyperbolic (torus) knots: {num_non_hyperbolic}")
    print(f"Number of hyperbolic knots: {num_hyperbolic}")
    print(f"The proportion of hyperbolic knots is the number of hyperbolic knots divided by the total number of knots.")
    print(f"Final equation: {num_hyperbolic} / {total_knots}")
    print(f"Result: {proportion}")

solve_knot_proportion()