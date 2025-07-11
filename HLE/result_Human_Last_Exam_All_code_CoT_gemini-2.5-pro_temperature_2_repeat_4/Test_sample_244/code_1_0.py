def calculate_hyperbolic_knot_proportion():
    """
    This function determines the proportion of hyperbolic knots among all knot types
    with a crossing number of exactly 7.
    """
    
    # Step 1: Identify all knot types with 7 crossings.
    # From standard knot tables (e.g., Rolfsen table or Knot Atlas), there are 7
    # distinct knot types with a crossing number of 7. They are denoted from 7_1 to 7_7.
    total_knots_count = 7
    
    # Step 2: Identify non-hyperbolic knots among them.
    # A prime knot is hyperbolic unless it's a torus knot. All 7-crossing knots are prime.
    # We check which of them are torus knots. The torus knot T(p,q) has a crossing
    # number of min(p(q-1), q(p-1)).
    # The only torus knot with 7 crossings is T(7,2), which is the knot 7_1.
    non_hyperbolic_count = 1  # This corresponds to the knot 7_1.
    
    # Step 3: Calculate the number of hyperbolic knots.
    # This is the total number of knots minus the non-hyperbolic ones.
    hyperbolic_count = total_knots_count - non_hyperbolic_count
    
    # Step 4: Calculate the proportion.
    proportion = hyperbolic_count / total_knots_count
    
    print(f"Total number of knot types with 7 crossings: {total_knots_count}")
    print(f"Number of non-hyperbolic (torus) knots: {non_hyperbolic_count}")
    print(f"Number of hyperbolic knots: {hyperbolic_count}")
    print(f"The proportion is calculated from the equation: {hyperbolic_count} / {total_knots_count}")
    print(f"As a decimal, the proportion is approximately: {proportion:.4f}")
    print(f"As a fraction, the proportion is: {hyperbolic_count}/{total_knots_count}")

calculate_hyperbolic_knot_proportion()
