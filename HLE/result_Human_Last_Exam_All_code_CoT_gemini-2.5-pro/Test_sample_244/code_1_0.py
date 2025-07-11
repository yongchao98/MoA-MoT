def calculate_hyperbolic_proportion():
    """
    Calculates the proportion of hyperbolic knots among those with 7 crossings.
    """
    # 1. Total number of knot types with 7 crossings.
    # According to the Rolfsen knot table, there are 7 prime knots with 7 crossings.
    total_knots = 7

    # 2. Number of non-hyperbolic knots with 7 crossings.
    # A prime knot is either a torus knot or a hyperbolic knot.
    # The only torus knot with 7 crossings is the 7_1 knot.
    non_hyperbolic_knots = 1

    # 3. Calculate the number of hyperbolic knots.
    hyperbolic_knots = total_knots - non_hyperbolic_knots

    # 4. Calculate the proportion.
    proportion = hyperbolic_knots / total_knots

    # Print the explanation and the final equation.
    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of non-hyperbolic (torus) knots: {non_hyperbolic_knots}")
    print(f"Number of hyperbolic knots: {hyperbolic_knots}")
    print("\nThe proportion of hyperbolic knots is calculated as:")
    print(f"Proportion = (Number of Hyperbolic Knots) / (Total Number of Knots)")
    print(f"Proportion = {hyperbolic_knots} / {total_knots}")
    print(f"\nThe decimal value of the proportion is: {proportion:.4f}")

calculate_hyperbolic_proportion()