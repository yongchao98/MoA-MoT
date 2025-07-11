def calculate_knot_proportion():
    """
    Calculates the proportion of 7-crossing knots that are hyperbolic.
    """
    # Step 1: Define the total number of knot types with 7 crossings.
    # The standard Rolfsen knot table lists 7 distinct prime knots for 7 crossings
    # when mirror images are considered the same type. These are denoted 7_1 through 7_7.
    total_knots = 7
    all_knots_list = ["7_1", "7_2", "7_3", "7_4", "7_5", "7_6", "7_7"]

    # Step 2: Identify the non-hyperbolic (i.e., torus) knots.
    # A prime knot is non-hyperbolic if and only if it is a torus knot.
    # For 7 crossings, the only torus knot is T(2,7), which is the knot 7_1.
    non_hyperbolic_knots = 1
    torus_knots_list = ["7_1"]

    # Step 3: Calculate the number of hyperbolic knots.
    # The hyperbolic knots are all knots that are not torus knots.
    hyperbolic_knots = total_knots - non_hyperbolic_knots

    # Step 4: Print the results and the calculation steps.
    print(f"Total number of distinct knot types with 7 crossings: {total_knots}")
    print(f"List of 7-crossing knots: {all_knots_list}")
    print("-" * 50)
    
    print(f"Number of non-hyperbolic (torus) knots: {non_hyperbolic_knots}")
    print(f"The non-hyperbolic knot is: {torus_knots_list[0]}")
    print("-" * 50)

    print("The number of hyperbolic knots is the total number of knots minus the non-hyperbolic knots.")
    print(f"Number of hyperbolic knots = {total_knots} - {non_hyperbolic_knots} = {hyperbolic_knots}")
    print("-" * 50)
    
    print("The proportion of hyperbolic knots is the number of hyperbolic knots divided by the total number of knots.")
    # The final equation with each number explicitly shown
    print(f"Proportion = {hyperbolic_knots} / {total_knots}")

    # Calculate and print the final proportion
    proportion = hyperbolic_knots / total_knots
    print(f"Final Proportion: {proportion:.4f}")

calculate_knot_proportion()