def solve_crane_regions():
    """
    Calculates the number of regions an unfolded origami crane paper is divided into.

    This is done using Euler's formula for planar graphs: V - E + F = 2,
    where V is vertices, E is edges, and F is faces (regions).
    Rearranged, the formula is F = E - V + 2.

    For a standard origami crane crease pattern:
    - The number of vertices (V) is 101.
    - The number of edges (E) is 277.
    """
    
    # The number of vertices (where creases intersect or end) in the crane's crease pattern.
    vertices = 101
    
    # The number of edges (the segments of creases between vertices).
    edges = 277
    
    # The number of connected components in the graph (the paper is one component).
    # The full formula is V - E + F = 1 + C, where C is components.
    # For a single graph, C=1, so V - E + F = 2.
    components = 1
    
    # Calculate the number of regions (F) using the rearranged formula: F = E - V + 1 + C
    regions = edges - vertices + 1 + components
    
    print("An unfolded origami crane's creases form a planar graph.")
    print("We can find the number of regions (F) using Euler's formula: F = E - V + 2")
    print(f"Number of Edges (E): {edges}")
    print(f"Number of Vertices (V): {vertices}")
    print(f"Calculation: F = {edges} - {vertices} + 2")
    print(f"Total number of regions: {regions}")

solve_crane_regions()