import networkx as nx

def create_super_knight_graph(n, m):
    """
    Constructs a (3,2) super-knight graph for an n x m board.
    """
    G = nx.Graph()
    # Add all squares as nodes
    nodes = [(r, c) for r in range(n) for c in range(m)]
    G.add_nodes_from(nodes)
    
    # Define the 8 possible moves for a (3,2) super-knight
    moves = [(3, 2), (3, -2), (-3, 2), (-3, -2),
             (2, 3), (2, -3), (-2, 3), (-2, -3)]

    # Add edges for valid moves
    for r1 in range(n):
        for c1 in range(m):
            for dr, dc in moves:
                r2, c2 = r1 + dr, c1 + dc
                # Check if the destination square is on the board
                if 0 <= r2 < n and 0 <= c2 < m:
                    G.add_edge((r1, c1), (r2, c2))
    return G

def find_max_planar_size():
    """
    Finds the supremum of nm for which the super-knight graph is planar.
    """
    max_planar_size = 0
    best_board_dims = None
    
    # We search n x m rectangles, with n<=m to avoid duplicates like 4x5 and 5x4.
    # The search range is chosen to be large enough to find the boundary.
    n_limit = 7
    m_limit = 12

    print("Checking board sizes (n x m) for planarity...")

    for n in range(4, n_limit):
        for m in range(n, m_limit):
            board_size = n * m
            g = create_super_knight_graph(n, m)
            is_planar, _ = nx.check_planarity(g)

            if is_planar:
                print(f"  {n}x{m} (size {board_size}): Planar")
                if board_size > max_planar_size:
                    max_planar_size = board_size
                    best_board_dims = (n, m)
            else:
                print(f"  {n}x{m} (size {board_size}): NON-PLANAR")
                # Since G(n,m) is non-planar, G(n, m+1), G(n, m+2), etc.
                # will also be non-planar, so we can stop searching for this n.
                break
    
    if best_board_dims:
        n, m = best_board_dims
        print("\n---")
        print(f"The largest planar board found is {n}x{m}.")
        # The final equation as requested
        print("The supremum of the value of the largest size is therefore:")
        print(f"{n} * {m} = {max_planar_size}")
    else:
        print("No planar graph found in the search range.")
        
# Run the search
find_max_planar_size()
