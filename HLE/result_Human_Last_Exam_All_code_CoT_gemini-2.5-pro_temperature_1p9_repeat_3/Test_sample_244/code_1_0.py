def solve_knot_proportion():
    """
    Calculates the proportion of hyperbolic knots among those with 7 crossings.

    This function relies on established data from knot theory:
    1. The number of distinct prime knots with a crossing number of 7 is 7.
       (These are denoted 7_1, 7_2, ..., 7_7).
    2. A prime knot is either a torus knot or a hyperbolic knot.
    3. Among knots with 7 crossings, only one is a torus knot: the 7_1 knot.
    """
    
    # Step 1: Define the total number of knot types with 7 crossings.
    # From knot tables, there are 7 such knots.
    total_knots = 7

    # Step 2: Define the number of non-hyperbolic (i.e., torus) knots.
    # For 7 crossings, only the 7_1 knot is a torus knot.
    non_hyperbolic_knots = 1

    # Step 3: Calculate the number of hyperbolic knots.
    # Hyperbolic knots are the total minus the non-hyperbolic ones.
    hyperbolic_knots = total_knots - non_hyperbolic_knots

    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of non-hyperbolic (torus) knots: {non_hyperbolic_knots}")
    print(f"Number of hyperbolic knots: {hyperbolic_knots}")
    print("-" * 30)
    print("The proportion of hyperbolic knots is the number of hyperbolic knots")
    print("divided by the total number of knots.")
    print("\nFinal Equation:")
    
    # Step 4: Output the numbers in the final equation.
    print(f"{hyperbolic_knots} / {total_knots}")

solve_knot_proportion()
