def solve_continuum_problem():
    """
    This program deduces the number of topologically distinct continua
    satisfying the given properties by following a logical argument.
    """

    print("Deduction Process:")
    print("1. Let X be the continuum, E be the set of its end points, and n = |E|.")
    print("   Property (1): 1 < n < infinity.")
    print("   Property (2): X has exactly two orbits, O1 and O2.")

    print("\n2. Since end points are a topological property, the set E must be an orbit.")
    print("   This means one orbit is E, and the other is X \\ E.")
    print("   - Orbit 1: The n end points.")
    print("   - Orbit 2: The infinite set of non-end-points.")

    print("\n3. All points in X \\ E are topologically equivalent. They must be cut-points.")
    print("   Let k be the number of components created by removing a point from X \\ E.")
    print("   For X \\ E to be a single orbit, k must be constant for all its points.")

    print("\n4. If k > 2, these are branch points. This would imply a complex structure with multiple point types (e.g., branch points and edge points), leading to more than one orbit in X \\ E. Thus, k must be 2.")
    print("   A space where every point is a cut-point of order 2 is homeomorphic to an open interval.")

    print("\n5. Therefore, X is a compactification of an open interval, with the n points of E added as the boundary.")
    print("   An open interval has two ends. A consistent compactification implies n=2.")
    print("   The resulting space is an arc (homeomorphic to [0, 1]).")

    print("\n6. Verification: An arc has 2 end points and 2 orbits ({endpoints}, {interior points}). It satisfies the conditions.")

    print("\nConclusion: Any such continuum must be homeomorphic to an arc. Therefore, there is only one such topological type.")

    # The final equation is simply the number of such continua.
    number_of_distinct_continua = 1
    
    print("\nFinal Answer Equation:")
    print(f"Number of topologically distinct continua = {number_of_distinct_continua}")

solve_continuum_problem()
<<<1>>>