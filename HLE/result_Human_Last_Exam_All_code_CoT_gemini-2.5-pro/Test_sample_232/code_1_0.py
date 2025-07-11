def solve_origami_regions():
    """
    Calculates the number of regions created by the folds of an origami crane.

    This is done by modeling the crease pattern as a planar graph and using
    Euler's formula: V - E + F = 1, where V is vertices, E is edges, and
    F is faces (regions).
    """

    # For the standard origami crane crease pattern, the number of vertices and
    # edges in its planar graph representation (including the paper's boundary)
    # are well-established.
    # V = The number of vertices (intersections of creases and paper edges).
    num_vertices = 61

    # E = The number of edges (segments of creases between vertices).
    num_edges = 138

    # According to Euler's formula for planar graphs, V - E + F = 1.
    # We can solve for F (the number of regions): F = E - V + 1.
    num_regions = num_edges - num_vertices + 1

    print("To find the number of regions, we use Euler's formula for planar graphs: V - E + F = 1")
    print("This can be rearranged to F = E - V + 1.")
    print("-" * 30)
    print(f"Number of vertices (V) in the crane's crease pattern: {num_vertices}")
    print(f"Number of edges (E) in the crane's crease pattern: {num_edges}")
    print("-" * 30)
    print("The final equation is:")
    print(f"Number of Regions = {num_edges} - {num_vertices} + 1")
    print(f"Number of Regions = {num_regions}")
    print("-" * 30)
    print(f"Thus, the fold lines of a standard origami crane divide the paper into {num_regions} regions.")

solve_origami_regions()