def solve_knot_proportion():
    """
    Calculates the proportion of 7-crossing knots that are hyperbolic.
    """
    # Step 1: Identify the total number of knot types with 7 crossings.
    # From knot theory tables (e.g., the Rolfsen knot table), there are 7
    # distinct prime knots with a crossing number of 7. These are denoted
    # 7_1, 7_2, 7_3, 7_4, 7_5, 7_6, and 7_7.
    total_knots = 7

    # Step 2: Identify the non-hyperbolic knots among them.
    # A prime knot is non-hyperbolic if and only if it is a torus knot.
    # We need to check which of the 7-crossing knots are torus knots.
    # The knot 7_1 is the torus knot T(7,2).
    # All other prime knots with 7 crossings (7_2 through 7_7) are not torus knots.
    non_hyperbolic_knots = 1 # This is the knot 7_1

    # Step 3: Calculate the number of hyperbolic knots.
    # The remaining knots are hyperbolic.
    hyperbolic_knots = total_knots - non_hyperbolic_knots

    # Step 4: Calculate and print the proportion.
    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of non-hyperbolic (torus) knots: {non_hyperbolic_knots}")
    print(f"Number of hyperbolic knots: {hyperbolic_knots}")
    print("\nThe proportion of these knots that are hyperbolic is the number of hyperbolic knots divided by the total number of knots.")
    print("The equation is:")
    print(f"{hyperbolic_knots} / {total_knots}")

solve_knot_proportion()