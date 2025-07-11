def solve_knot_proportion():
    """
    Calculates the proportion of hyperbolic knots among all knot types
    with exactly 7 crossings.
    """
    # Step 1 & 2: List the 7 prime knots with 7 crossings.
    prime_knots_7_crossings = ["7_1", "7_2", "7_3", "7_4", "7_5", "7_6", "7_7"]
    num_prime_knots = len(prime_knots_7_crossings)

    # Step 3: List the composite knots with 7 crossings.
    # The only combination is 3_1 + 4_1.
    composite_knots_7_crossings = ["3_1#4_1"]
    num_composite_knots = len(composite_knots_7_crossings)

    # Step 4: Find the total number of knots.
    total_knots_list = prime_knots_7_crossings + composite_knots_7_crossings
    total_knots = len(total_knots_list)

    print(f"There are {num_prime_knots} prime knots and {num_composite_knots} composite knot(s) with 7 crossings.")
    print(f"Total number of knot types with 7 crossings: {num_prime_knots} + {num_composite_knots} = {total_knots}")
    print("The list of all 7-crossing knots is:", ", ".join(total_knots_list))
    print("-" * 20)

    # Step 5: Identify which of these are hyperbolic.
    # A knot is hyperbolic if it is not a torus knot or a satellite knot.
    # 7_1 is the (7,2)-torus knot.
    # 3_1#4_1 is a composite knot, which is a type of satellite knot.
    # The remaining knots (7_2 to 7_7) are hyperbolic.
    torus_knots = ["7_1"]
    satellite_knots = ["3_1#4_1"]
    non_hyperbolic_knots = torus_knots + satellite_knots

    hyperbolic_knots = [k for k in total_knots_list if k not in non_hyperbolic_knots]
    num_hyperbolic = len(hyperbolic_knots)
    
    print(f"Among these, the non-hyperbolic knots are: {', '.join(non_hyperbolic_knots)}")
    print(f"The number of hyperbolic knots is: {num_hyperbolic}")
    print("The list of hyperbolic knots is:", ", ".join(hyperbolic_knots))
    print("-" * 20)

    # Step 6: Calculate and print the proportion.
    if total_knots > 0:
        proportion = num_hyperbolic / total_knots
        print("The proportion of hyperbolic knots is the number of hyperbolic knots divided by the total number of knots.")
        print("Final Calculation:")
        print(f"{num_hyperbolic} / {total_knots} = {proportion}")
    else:
        print("There are no knots to consider.")

solve_knot_proportion()
<<<0.75>>>