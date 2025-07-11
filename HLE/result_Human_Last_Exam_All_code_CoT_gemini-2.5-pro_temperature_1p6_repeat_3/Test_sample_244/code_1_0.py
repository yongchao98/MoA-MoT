def calculate_hyperbolic_knot_proportion():
    """
    Calculates the proportion of hyperbolic knots among those with 7 crossings.

    According to knot theory, there are 7 distinct types of knots with exactly 7 crossings.
    Of these, only one, the 7_1 knot (which is the (7,2)-torus knot), is not hyperbolic.
    The remaining 6 knots (7_2 through 7_7) are hyperbolic.
    """
    total_knots = 7
    non_hyperbolic_knots = 1 # This is the 7_1 knot (a torus knot)
    
    hyperbolic_knots = total_knots - non_hyperbolic_knots
    
    # Calculate the proportion
    proportion = hyperbolic_knots / total_knots
    
    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of hyperbolic knots with 7 crossings: {hyperbolic_knots}")
    print(f"The proportion is calculated as the ratio of hyperbolic knots to the total number of knots.")
    print(f"Equation: {hyperbolic_knots} / {total_knots} = {proportion}")

if __name__ == "__main__":
    calculate_hyperbolic_knot_proportion()
