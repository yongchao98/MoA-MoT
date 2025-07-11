def calculate_hyperbolic_knot_proportion():
    """
    This function calculates the proportion of hyperbolic knots among all knot
    types with exactly 7 crossings.

    According to knot theory, there are 7 distinct prime knots with 7 crossings:
    7_1, 7_2, 7_3, 7_4, 7_5, 7_6, and 7_7.

    Of these, only one is not hyperbolic:
    - 7_1 is the (7,2)-torus knot. Torus knots are not hyperbolic.

    The other 6 knots are all hyperbolic.
    """
    
    # Total number of distinct knot types with 7 crossings
    total_knots = 7
    
    # Number of hyperbolic knots among them
    num_hyperbolic_knots = 6
    
    # Calculate the proportion
    proportion = num_hyperbolic_knots / total_knots
    
    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of hyperbolic knots with 7 crossings: {num_hyperbolic_knots}")
    print("\nThe proportion is calculated as:")
    print(f"{num_hyperbolic_knots} / {total_knots} = {proportion}")

calculate_hyperbolic_knot_proportion()