import spherogram

def solve_knot_proportion():
    """
    Calculates the proportion of hyperbolic knots among all knots
    with a specific crossing number.
    """
    crossing_number = 7
    
    # Spherogram's census contains prime knots, and by convention,
    # a knot and its mirror are considered the same type.
    try:
        knots_at_crossing = spherogram.knot_census(crossing_number)
    except (IndexError, KeyError):
        print(f"No knot data available for {crossing_number} crossings in the census.")
        return

    total_knots = 0
    hyperbolic_knots = 0

    for knot in knots_at_crossing:
        total_knots += 1
        # A knot is hyperbolic if its complement admits a complete
        # hyperbolic metric of finite volume. For small crossing numbers,
        # all prime knots are hyperbolic.
        if knot.is_hyperbolic():
            hyperbolic_knots += 1
    
    print(f"Considering all knot types with exactly {crossing_number} crossings:")
    print(f"Total number of distinct knot types = {total_knots}")
    print(f"Number of hyperbolic knot types = {hyperbolic_knots}")
    
    if total_knots > 0:
        proportion = hyperbolic_knots / total_knots
        # The final equation is printed here
        print(f"The proportion of hyperbolic knots is {hyperbolic_knots} / {total_knots} = {proportion}")
    else:
        print("Cannot calculate proportion as there are no knots.")

solve_knot_proportion()
