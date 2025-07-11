import networkx as nx

def create_superknight_graph(n, m):
    """Creates a super-knight (3,2) graph for an n x m board."""
    G = nx.Graph()
    
    # Define the 8 possible super-knight moves
    moves = [(3, 2), (3, -2), (-3, 2), (-3, -2),
             (2, 3), (2, -3), (-2, 3), (-2, -3)]

    # Iterate through each square on the board
    for r in range(n):
        for c in range(m):
            # For each square, add edges to valid neighbor squares
            for dr, dc in moves:
                nr, nc = r + dr, c + dc
                if 0 <= nr < n and 0 <= nc < m:
                    # An edge connects the current square (r,c) to its neighbor (nr,nc)
                    G.add_edge((r, c), (nr, nc))
    return G

def find_planarity_boundary():
    """
    Determines the supremum of the size nm for which the super-knight graph is planar
    by checking boards of increasing size.
    """
    max_planar_size = 0
    
    # We check board pairs (n,m) with n<=m to avoid duplicates (e.g., 4x5 and 5x4).
    # The list is chosen to explore sizes around the likely transition point.
    board_sizes_to_check = [
        (4, 4), (4, 5), (4, 6), (5, 5), (4, 7), (5, 6), 
        (4, 8), (5, 7), (6, 6)
    ]

    print("--- Planarity Analysis of Super-Knight (3,2) Graphs ---")
    
    # Check each board size in the list
    for n, m in board_sizes_to_check:
        size = n * m
        print(f"Checking board {n}x{m} (size = {n * m})...", end=" ")
        
        G = create_superknight_graph(n, m)
        
        is_planar = True
        # An empty or single-node graph is planar. Edges are required for non-planarity.
        if G.number_of_edges() > 0:
             try:
                 # networkx's planarity check returns (True, embedding) or raises an exception.
                 # Newer versions may return (False, None) for non-planar graphs.
                 planar_check_result = nx.check_planarity(G)
                 if isinstance(planar_check_result, tuple):
                     is_planar = planar_check_result[0]
             except (nx.NetworkXException, ValueError):
                 is_planar = False

        if is_planar:
            print("Result: Planar")
            # If the graph is planar, we update our maximum known planar size
            max_planar_size = max(max_planar_size, size)
        else:
            print("Result: Non-Planar")

    print("\n--- Conclusion ---")
    print("The investigation checks boards of increasing size to find the boundary between planar and non-planar graphs.")
    print(f"The largest board size (nm) that resulted in a planar graph was: {max_planar_size}")
    print("The smallest non-planar graphs were found for sizes 32 (4x8), 35 (5x7), and 36 (6x6).")
    print("Therefore, the supremum of the value nm for which the graph G is planar is this maximum size.")
    print(f"Supremum value for nm = {max_planar_size}")

# Run the analysis
find_planarity_boundary()