import math

def solve_knot_proportion():
    """
    This function calculates the proportion of hyperbolic knots among all knot types
    with exactly 7 crossings.

    According to knot theory and the Alexander-Briggs notation, there are 7 distinct
    knot types (up to mirroring) with 7 crossings. These are denoted 7_1 through 7_7.

    A prime knot is either a torus knot or a hyperbolic knot. We need to identify
    which of the 7-crossing knots are hyperbolic.
    - The knot 7_7 is known to be the (7,2)-torus knot. Therefore, it is not hyperbolic.
    - The other six knots (7_1, 7_2, 7_3, 7_4, 7_5, and 7_6) are not torus knots,
      and since they are prime, they are all hyperbolic.

    Thus, we can calculate the proportion.
    """

    # Total number of knot types with 7 crossings
    total_knots = 7

    # Number of hyperbolic knot types with 7 crossings
    # (7_1, 7_2, 7_3, 7_4, 7_5, 7_6)
    hyperbolic_knots = 6

    # Calculate the proportion
    proportion = hyperbolic_knots / total_knots

    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of hyperbolic knot types among them: {hyperbolic_knots}")
    print("\nThe proportion of hyperbolic knots is calculated as:")
    # We output each number in the final equation.
    print(f"{hyperbolic_knots} / {total_knots} = {proportion}")

solve_knot_proportion()