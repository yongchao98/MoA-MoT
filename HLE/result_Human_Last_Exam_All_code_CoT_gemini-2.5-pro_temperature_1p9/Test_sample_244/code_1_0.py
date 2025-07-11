import spherogram

def solve_knot_proportion():
    """
    Calculates the proportion of hyperbolic knots among those with 7 crossings.
    """
    # Using spherogram's built-in knot catalog to get all knots with 7 crossings.
    # The standard catalog considers a knot and its mirror as the same type.
    knots_7_crossings = spherogram.knots.identify('K7')

    # Get the total count of these knots.
    total_knots = len(knots_7_crossings)

    # A prime knot is hyperbolic if and only if it is not a torus knot.
    # All knots with 7 crossings are prime.
    # We identify the hyperbolic ones by checking which are not torus knots.
    hyperbolic_knots_count = 0
    for knot in knots_7_crossings:
        if not knot.is_torus():
            hyperbolic_knots_count += 1

    # The proportion is the number of hyperbolic knots divided by the total number of knots.
    proportion = hyperbolic_knots_count / total_knots

    print(f"There are {total_knots} distinct knot types with 7 crossings.")
    print(f"Out of these, {hyperbolic_knots_count} are hyperbolic.")
    print("The proportion of hyperbolic knots is calculated as follows:")
    print(f"{hyperbolic_knots_count} / {total_knots} = {proportion}")

solve_knot_proportion()