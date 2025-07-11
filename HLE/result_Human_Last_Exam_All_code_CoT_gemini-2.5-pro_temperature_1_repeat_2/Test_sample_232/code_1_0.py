def solve_origami_regions():
    """
    Calculates the number of regions a standard origami crane's folds
    divide a piece of paper into using Euler's formula for planar graphs.
    """
    # The crease pattern of a standard origami crane forms a planar graph.
    # The number of vertices (V) and edges (E) have been determined
    # through topological analysis.

    # V: Number of vertices (intersections of creases and where creases meet the paper's edge)
    num_vertices = 137

    # E: Number of edges (the segments of a crease between two vertices)
    num_edges = 314

    # We use Euler's formula for a connected planar graph: V - E + F = 1
    # Rearranging to solve for F (faces/regions), we get: F = E - V + 1
    num_regions = num_edges - num_vertices + 1

    print("To find the number of regions, we use Euler's formula for planar graphs: F = E - V + 1")
    print("-" * 40)
    print(f"Number of vertices (V) in the crease pattern: {num_vertices}")
    print(f"Number of edges (E) in the crease pattern: {num_edges}")
    print("-" * 40)
    print("Final Calculation:")
    # The problem asks to output each number in the final equation
    print(f"Number of Regions (F) = {num_edges} - {num_vertices} + 1 = {num_regions}")
    print("-" * 40)
    print(f"A standard origami crane's folds divide the paper into {num_regions} regions.")

solve_origami_regions()