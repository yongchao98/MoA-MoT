import math

def solve_knot_proportion():
    """
    Calculates the proportion of hyperbolic knots among those with exactly 7 crossings.
    """
    # Step 1: Identify the total number of knot types with 7 crossings.
    # From the Rolfsen knot table, there are 7 distinct prime knots with 7 crossings.
    # (7_1, 7_2, 7_3, 7_4, 7_5, 7_6, 7_7)
    total_knots = 7

    # Step 2: Identify the number of non-hyperbolic knots with 7 crossings.
    # A prime knot is non-hyperbolic if it is a torus or a satellite knot.
    # Based on knot theory classification:
    # - The torus knots with fewer than 8 crossings are 3_1 and 5_1. There are none with 7 crossings.
    # - The first satellite knot (8_19) has 8 crossings. There are none with 7 crossings.
    num_torus_knots = 0
    num_satellite_knots = 0
    num_non_hyperbolic_knots = num_torus_knots + num_satellite_knots

    # Step 3: Calculate the number of hyperbolic knots.
    num_hyperbolic_knots = total_knots - num_non_hyperbolic_knots

    # Step 4: Calculate the proportion.
    if total_knots == 0:
        proportion = 0.0
    else:
        proportion = num_hyperbolic_knots / total_knots

    # Step 5: Print the results and the final equation.
    print(f"Total number of knot types with 7 crossings: {total_knots}")
    print(f"Number of non-hyperbolic knots (torus or satellite): {num_non_hyperbolic_knots}")
    print(f"Number of hyperbolic knots: {num_hyperbolic_knots}")
    print("The proportion of hyperbolic knots is calculated as:")
    print(f"Proportion = (Number of Hyperbolic Knots) / (Total Number of Knots)")
    print(f"Proportion = {num_hyperbolic_knots} / {total_knots} = {proportion}")

solve_knot_proportion()