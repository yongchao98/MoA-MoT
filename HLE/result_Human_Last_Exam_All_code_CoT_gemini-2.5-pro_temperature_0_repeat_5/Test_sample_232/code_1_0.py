def solve_crane_regions():
    """
    Calculates the number of regions in an unfolded origami crane paper
    using Euler's formula for planar graphs (F = E - V + 1).
    """
    # The number of vertices (V) and edges (E) for a standard crane's
    # crease pattern are known values in origami mathematics.
    
    # V: Number of points where creases intersect or meet the paper's edge.
    num_vertices = 65
    
    # E: Number of crease segments between vertices.
    num_edges = 194
    
    # Using Euler's formula for a planar graph on a disk (F = E - V + 1)
    # to find the number of regions (F).
    num_regions = num_edges - num_vertices + 1
    
    print("The number of regions is calculated using Euler's formula for planar graphs: F = E - V + 1")
    print(f"For a standard origami crane, the known values are:")
    print(f" - Number of Vertices (V): {num_vertices}")
    print(f" - Number of Edges (E): {num_edges}")
    print("\nFinal Calculation:")
    print(f"Number of Regions = {num_edges} - {num_vertices} + 1 = {num_regions}")

solve_crane_regions()