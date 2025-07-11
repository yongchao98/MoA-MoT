def calculate_hyperbolic_knot_proportion():
    """
    Calculates the proportion of hyperbolic knots among all knot types
    with a crossing number of 7.
    """
    # 1. Identify all prime and composite knots with 7 crossings.
    # There are 7 prime knots with 7 crossings.
    prime_knots_7_crossings = ['7_1', '7_2', '7_3', '7_4', '7_5', '7_6', '7_7']
    num_prime_knots = len(prime_knots_7_crossings)

    # The only composite knot with 7 crossings is the sum of the 3_1 and 4_1 knots.
    composite_knots_7_crossings = ['3_1 # 4_1']
    num_composite_knots = len(composite_knots_7_crossings)

    # 2. Calculate the total number of knots.
    total_knots = num_prime_knots + num_composite_knots

    # 3. Identify the non-hyperbolic knots.
    # These are torus knots and satellite knots (which include all composite knots).
    # From the list of 7-crossing knots:
    # - 7_1 is the torus knot T(7,2).
    # - 3_1 # 4_1 is a composite knot.
    non_hyperbolic_knots = ['7_1', '3_1 # 4_1']
    num_non_hyperbolic = len(non_hyperbolic_knots)

    # 4. Calculate the number of hyperbolic knots.
    num_hyperbolic = total_knots - num_non_hyperbolic

    # 5. Calculate and print the proportion.
    proportion = num_hyperbolic / total_knots

    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of hyperbolic knots: {num_hyperbolic}")
    print(f"Number of non-hyperbolic knots: {num_non_hyperbolic}")
    print("\nThe proportion of hyperbolic knots is calculated as:")
    # Final equation with numbers
    print(f"Proportion = {num_hyperbolic} / {total_knots}")
    print(f"Result: {proportion}")

calculate_hyperbolic_knot_proportion()