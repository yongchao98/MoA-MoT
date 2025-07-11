def solve_knot_proportion():
    """
    Calculates the proportion of 7-crossing knots that are hyperbolic.
    
    This is based on established knot theory classifications.
    """
    
    # Step 1: Define the total number of distinct knot types with 7 crossings.
    # According to knot tables, these are 7_1 through 7_7.
    total_knots = 7
    
    # Step 2: Identify the number of non-hyperbolic knots.
    # For prime knots, non-hyperbolic knots are torus knots.
    # Among the 7-crossing knots, only 7_7 (the T(7,2) torus knot) is non-hyperbolic.
    non_hyperbolic_knots = 1
    
    # Step 3: Calculate the number of hyperbolic knots.
    hyperbolic_knots = total_knots - non_hyperbolic_knots
    
    # Step 4: Display the result as a descriptive equation.
    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of hyperbolic knots: {hyperbolic_knots}")
    print("\nThe proportion of hyperbolic knots is calculated as:")
    print(f"{hyperbolic_knots} / {total_knots}")

solve_knot_proportion()