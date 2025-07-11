import spherogram

def solve_knot_proportion():
    """
    Calculates the proportion of 7-crossing knots that are hyperbolic.
    
    A prime knot is hyperbolic if and only if it is not a torus knot.
    The standard knot tables K7_1, ..., K7_7 list the 7 distinct prime
    knots with 7 crossings.
    """
    # List of all knot types with 7 crossings
    knot_identifiers = [f'7_{i}' for i in range(1, 8)]
    
    total_knots = 0
    hyperbolic_knots = 0
    
    # Check each knot
    for ident in knot_identifiers:
        knot = spherogram.Knot(ident)
        total_knots += 1
        
        # A prime knot is hyperbolic if and only if it is not a torus knot.
        if not knot.is_torus():
            hyperbolic_knots += 1
            
    if total_knots == 0:
        proportion = 0
    else:
        proportion = hyperbolic_knots / total_knots

    print(f"There are {total_knots} distinct knot types with 7 crossings.")
    print(f"Among these, {hyperbolic_knots} are hyperbolic.")
    print("The proportion of hyperbolic knots is calculated as follows:")
    print(f"{hyperbolic_knots} / {total_knots} = {proportion}")

solve_knot_proportion()
<<<0.8571428571428571>>>