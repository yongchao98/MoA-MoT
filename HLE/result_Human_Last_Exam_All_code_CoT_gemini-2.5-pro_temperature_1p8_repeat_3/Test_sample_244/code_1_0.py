def solve_knot_proportion():
    """
    Calculates the proportion of hyperbolic knots among those with 7 crossings.

    This is based on established knot theory classifications.
    """

    # According to knot theory tables, there are 7 distinct knot types with a minimal crossing number of 7.
    total_knots = 7

    # Of these 7 knots, only one (the torus knot 7_1) is not hyperbolic.
    # The remaining 6 knots (7_2, 7_3, 7_4, 7_5, 7_6, 7_7) are hyperbolic.
    num_hyperbolic = 6
    
    # Calculate the proportion
    proportion = num_hyperbolic / total_knots

    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of hyperbolic knots among them: {num_hyperbolic}")
    print("\nThe proportion is the number of hyperbolic knots divided by the total number of knots.")
    print(f"Final Equation: {num_hyperbolic} / {total_knots} = {proportion}")

solve_knot_proportion()