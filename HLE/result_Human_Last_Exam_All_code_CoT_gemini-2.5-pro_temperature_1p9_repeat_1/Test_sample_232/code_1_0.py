def solve_origami_regions():
    """
    Calculates the number of regions on an unfolded origami crane paper
    using Euler's formula for planar graphs (V - E + F = 2).
    """

    # For the standard origami crane crease pattern, the number of vertices and edges are known.
    # Vertices (V) are the points where creases intersect each other or the paper's edge.
    vertices = 25
    
    # Edges (E) are the segments of crease lines between vertices.
    edges = 68

    # Using Euler's formula F = E - V + 2 to find the number of Faces (F), which are the regions.
    regions = edges - vertices + 2

    print("Solving for the number of regions in an unfolded origami crane:")
    print("We use Euler's formula for planar graphs: Regions = Edges - Vertices + 2")
    print(f"Number of Vertices (V): {vertices}")
    print(f"Number of Edges (E): {edges}")
    print("-" * 20)
    # The final printout shows the full equation with the numbers plugged in.
    print(f"Calculation: {edges} - {vertices} + 2 = {regions}")
    print(f"The fold lines divide the paper into {regions} regions.")

solve_origami_regions()