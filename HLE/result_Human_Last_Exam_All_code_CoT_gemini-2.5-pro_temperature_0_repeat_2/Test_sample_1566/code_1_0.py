def solve_continuum_problem():
    """
    This function solves the topological problem by analyzing the given properties
    and identifying the continua that satisfy them.
    """

    # Let N be the number of topologically distinct continua satisfying the properties.

    # Property (1): The continuum X has n end points, where 1 < n < infinity.
    # Property (2): X has exactly 2 orbits under auto-homeomorphisms, O_1 and O_2.

    # The set of end points, E, must be one of the orbits. Let E = O_1.
    # This means |O_1| = n, and O_2 is the set of all non-end points.
    # The homogeneity of the orbit O_2 (the non-end points) strongly restricts the
    # structure of X. A detailed analysis of the local topology shows that the
    # non-end points must all be locally like the interior of an arc.
    # This forces O_2 to be topologically an open interval.

    # The space X must be a compactification of this open interval. The points
    # added to form the compact space must be the end points, O_1.
    # The standard compactification of an open interval is the closed interval,
    # which adds two points at the ends. This results in a space that is
    # topologically equivalent to the closed interval [0, 1].

    # Let's verify this candidate: X = [0, 1].
    # 1. End points: The points 0 and 1 are the end points. The number of end points is 2.
    #    This satisfies the condition "more than one and finitely many".
    # 2. Orbits: The homeomorphisms of [0, 1] create exactly two orbits:
    #    - The set of end points: {0, 1}
    #    - The set of interior points: (0, 1)
    #    This satisfies the condition of having exactly two orbits.

    # The analysis suggests that the closed interval is the only topological space
    # that meets these stringent requirements.

    # Therefore, there is only one such topologically distinct continuum.
    number_of_continua = 1

    # The problem asks for the number of such continua.
    # The final equation is simply: Total = 1
    print("The analysis leads to a unique topological space (the closed interval).")
    print(f"The number of topologically distinct continua is: {number_of_continua}")

solve_continuum_problem()