def calculate_hyperbolic_knot_proportion():
    """
    Calculates the proportion of 7-crossing knots that are hyperbolic.

    A knot is hyperbolic if it is a prime knot and not a torus knot.
    We use standard knot theory classifications for knots with 7 crossings.
    """

    # According to knot theory, there are 7 distinct knots with 7 crossings.
    # We represent them here and mark whether they are hyperbolic (True) or not (False).
    # All prime knots that are not torus knots are hyperbolic.
    knots_7_crossings = {
        "7_1": True,  # Prime, not a torus knot -> Hyperbolic
        "7_2": True,  # Prime, not a torus knot -> Hyperbolic
        "7_3": True,  # Prime, not a torus knot -> Hyperbolic
        "7_4": True,  # Prime, not a torus knot -> Hyperbolic
        "7_5": True,  # Prime, not a torus knot -> Hyperbolic
        "7_6": True,  # Prime, not a torus knot -> Hyperbolic
        "7_7": False, # This is the (7,2)-torus knot -> Not Hyperbolic
    }

    # 1. Count the total number of knots
    total_knots = len(knots_7_crossings)

    # 2. Count the number of hyperbolic knots
    # This is done by summing up the 'True' values in our dictionary
    num_hyperbolic_knots = sum(knots_7_crossings.values())

    # 3. Calculate the proportion
    proportion = num_hyperbolic_knots / total_knots

    # Print the results and the final equation
    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of hyperbolic knots among them: {num_hyperbolic_knots}")
    print("\nThe proportion is calculated as:")
    print(f"Proportion = Number of Hyperbolic Knots / Total Number of Knots")
    print(f"Equation: {num_hyperbolic_knots} / {total_knots}")
    print(f"\nThe result is: {proportion}")


calculate_hyperbolic_knot_proportion()
<<<0.8571428571428571>>>