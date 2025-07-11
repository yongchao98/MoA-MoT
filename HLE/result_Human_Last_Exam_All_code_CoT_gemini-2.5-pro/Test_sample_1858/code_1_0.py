def solve_polygon_components():
    """
    Calculates the number of connected components of the space of
    non-self-intersecting 6-sided polygons in R^3.

    The solution is based on principles from knot theory.
    """

    # Step 1: Decompose the space by knot type.
    # The space of non-self-intersecting polygons is partitioned by the knot type of the polygon.
    # A continuous deformation (an isotopy) preserves the knot type.
    # Thus, the total number of connected components is the sum of the number of components for each possible knot type.

    # Step 2: Identify possible knot types for a 6-sided polygon.
    # This is determined by the "stick number" of a knot, which is the minimum
    # number of segments needed to construct it.
    # Stick numbers for the simplest knots:
    # - Unknot (0_1): 3
    # - Trefoil knot (3_1): 6
    # - Figure-eight knot (4_1): 7
    # A 6-sided polygon can only form knots with a stick number of 6 or less.
    # Therefore, a 6-sided polygon can only be an unknot or a trefoil knot.

    # Step 3: Count the connected components for each knot type.

    # For the unknot:
    # The space of n-sided unknotted polygons is known to be path-connected.
    # This means all 6-sided unknots can be deformed into one another.
    num_components_unknot = 1

    # For the trefoil knot:
    # The trefoil knot is chiral; it exists in right-handed and left-handed forms
    # which are not isotopic to each other. One cannot be deformed into the other
    # without self-intersection.
    # The space of 6-sided trefoil knots consists of two connected components.
    num_components_trefoil = 2

    # Step 4: Calculate the total number of connected components.
    total_components = num_components_unknot + num_components_trefoil

    print("The total number of connected components is the sum of components for each possible knot type.")
    print(f"Number of components for 6-sided unknots: {num_components_unknot}")
    print(f"Number of components for 6-sided trefoil knots (right- and left-handed): {num_components_trefoil}")
    print(f"Total number of connected components = {num_components_unknot} + {num_components_trefoil} = {total_components}")

solve_polygon_components()