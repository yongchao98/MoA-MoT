import math

def calculate_hyperbolic_knot_proportion():
    """
    Calculates the proportion of hyperbolic knots among those with 7 crossings.
    """
    # Step 1: Find the total number of knot types with 7 crossings.
    # From standard knot tables (e.g., Rolfsen table), the knots with 7 crossings
    # are 7_1, 7_2, 7_3, 7_4, 7_5, 7_6, and 7_7.
    total_knots = 7
    knot_types = ["7_1", "7_2", "7_3", "7_4", "7_5", "7_6", "7_7"]

    # Step 2 & 3: Identify the non-hyperbolic (torus) knots.
    # A prime knot is either a torus knot or a hyperbolic knot. All knots in the
    # list are prime. The only torus knot among them is the 7_1 knot, which is
    # the T(2,7) torus knot.
    non_hyperbolic_knots_list = ["7_1"]
    num_non_hyperbolic_knots = len(non_hyperbolic_knots_list)

    # The rest are hyperbolic.
    num_hyperbolic_knots = total_knots - num_non_hyperbolic_knots

    # Step 4: Calculate the proportion.
    proportion = num_hyperbolic_knots / total_knots

    print(f"Total number of distinct knot types with 7 crossings: {total_knots}")
    print(f"These knots are: {', '.join(knot_types)}")
    print(f"Among these, the non-hyperbolic (torus) knot is: {', '.join(non_hyperbolic_knots_list)}")
    print(f"Therefore, the number of hyperbolic knots is: {num_hyperbolic_knots}")
    print("\nThe proportion of hyperbolic knots is calculated as:")
    print(f"Equation: {num_hyperbolic_knots} / {total_knots}")
    print(f"Decimal value: {proportion:.4f}")

calculate_hyperbolic_knot_proportion()