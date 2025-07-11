def solve_knot_proportion():
    """
    Calculates the proportion of 7-crossing knots that are hyperbolic.
    """
    # Step 1: Total number of knot types with 7 crossings.
    # From knot theory tables (e.g., Rolfsen's knot table), there are 7
    # distinct knot types with 7 crossings, denoted 7_1 through 7_7.
    total_knots = 7

    # Step 2: Identify non-hyperbolic (torus) knots.
    # A prime knot is non-hyperbolic if and only if it is a torus knot.
    # All 7-crossing knots are prime.
    # The crossing number of a (p,q)-torus knot is given by min(p(q-1), q(p-1)).
    # For a crossing number of 7, the only solution for coprime integers p, q > 1
    # is the (2,7)-torus knot. This knot is conventionally denoted as 7_1.
    num_torus_knots = 1

    # Step 3: Calculate the number of hyperbolic knots.
    # The remaining knots are hyperbolic.
    num_hyperbolic_knots = total_knots - num_torus_knots

    # Step 4: Calculate the proportion.
    proportion = num_hyperbolic_knots / total_knots

    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of non-hyperbolic (torus) knots: {num_torus_knots}")
    print(f"Number of hyperbolic knots: {num_hyperbolic_knots}")
    print("\nThe proportion of hyperbolic knots is calculated as:")
    print(f"{num_hyperbolic_knots} / {total_knots} = {proportion:.5f}")

solve_knot_proportion()