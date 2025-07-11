def solve_knot_components():
    """
    Calculates the number of connected components of the space of
    non-self-intersecting 6-sided polygons in R^3.

    This is equivalent to finding the number of knot types with a stick number
    of 6 or less.
    """

    # The maximum number of sticks allowed by the polygon's sides.
    max_sticks = 6

    # --- Analysis of Knot Types ---

    # 1. The Unknot (0_1)
    # The stick number of the unknot is 3.
    s_unknot = 3
    # Check if it can be formed with the given number of sticks.
    if s_unknot <= max_sticks:
        num_unknot_components = 1
    else:
        num_unknot_components = 0
    
    # 2. The Trefoil Knot (3_1)
    # The stick number of the trefoil knot is 6.
    s_trefoil = 6
    # The trefoil is chiral, so it has two distinct types (left- and right-handed).
    if s_trefoil <= max_sticks:
        num_trefoil_components = 2
    else:
        num_trefoil_components = 0

    # 3. The Figure-Eight Knot (4_1)
    # The stick number of the figure-eight knot is 7.
    s_figure_eight = 7
    if s_figure_eight <= max_sticks:
        num_figure_eight_components = 1
    else:
        num_figure_eight_components = 0
    
    # Other knots have stick numbers greater than 6.
    
    # --- Final Calculation ---
    total_components = num_unknot_components + num_trefoil_components + num_figure_eight_components
    
    print("The total number of connected components is the sum of possible knot types.")
    print(f"1. Unknot: Stick number is {s_unknot}. Possible with <= {max_sticks} sticks. Components: {num_unknot_components}")
    print(f"2. Trefoil Knots: Stick number is {s_trefoil}. Possible with <= {max_sticks} sticks. Components: {num_trefoil_components}")
    print(f"3. Figure-Eight Knot: Stick number is {s_figure_eight}. Not possible with <= {max_sticks} sticks. Components: {num_figure_eight_components}")
    print("\nFinal Equation:")
    print(f"{num_unknot_components} + {num_trefoil_components} = {total_components}")

solve_knot_components()