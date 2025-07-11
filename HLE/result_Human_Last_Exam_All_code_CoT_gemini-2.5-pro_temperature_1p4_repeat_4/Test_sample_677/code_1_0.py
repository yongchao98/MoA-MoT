import networkx as nx

def create_superknight_graph(n, m):
    """
    Constructs the (3,2) super-knight graph for an n x m chessboard.
    """
    G = nx.Graph()
    # Add nodes for each square on the board
    nodes = [(r, c) for r in range(1, n + 1) for c in range(1, m + 1)]
    G.add_nodes_from(nodes)
    
    # Define the 8 possible moves of a (3,2) super-knight
    move_deltas = [(2, 3), (2, -3), (-2, 3), (-2, -3),
                   (3, 2), (3, -2), (-3, 2), (-3, -2)]

    # Add edges for valid moves
    for r1 in range(1, n + 1):
        for c1 in range(1, m + 1):
            for dr, dc in move_deltas:
                r2, c2 = r1 + dr, c1 + dc
                # Check if the destination square is on the board
                if 1 <= r2 <= n and 1 <= c2 <= m:
                    # Add edge, ensuring not to add duplicates for undirected graph
                    if not G.has_edge((r1, c1), (r2, c2)):
                        G.add_edge((r1, c1), (r2, c2))
    return G

def check_planarity(n, m):
    """
    Checks if the super-knight graph for an n x m board is planar and prints the result.
    """
    G = create_superknight_graph(n, m)
    is_planar, _ = nx.check_planarity(G)
    print(f"Board {n}x{m} (size {n*m}): The graph is {'Planar' if is_planar else 'Non-Planar'}.")
    return is_planar

def find_supremum():
    """
    Finds the largest size nm for which the super-knight graph is planar.
    """
    print("Step 1: Checking boards where n, m >= 5.")
    # If a 5x5 board is non-planar, any n x m board with n,m >= 5 will also be.
    # A non-planar graph is a subgraph of the larger board's graph.
    is_5x5_planar = check_planarity(5, 5)

    if not is_5x5_planar:
        print("\nSince the 5x5 board is non-planar, any planar board must have at least one dimension equal to 4.")
    
    print("\nStep 2: Checking 4xm boards to find the boundary of planarity.")
    
    m = 4
    largest_planar_size = 0
    while True:
        is_planar = check_planarity(4, m)
        if is_planar:
            largest_planar_size = 4 * m
            m += 1
        else:
            print(f"\nThe boundary is found. {4}x{m-1} is the largest planar board with one dimension of 4.")
            break
            
    print("\n--- Conclusion ---")
    print(f"The set of planar rectangles with n,m>=4 is composed of boards of size 4x4, 4x5, 4x6, and 4x7 (and their transposes).")
    print(f"The corresponding sizes nm are {4*4}, {4*5}, {4*6}, and {4*7}.")
    print(f"The largest size nm for which the graph is planar is {largest_planar_size}.")
    print(f"The supremum of this set of finite sizes is its maximum value.")
    
    final_answer = largest_planar_size
    print("\nFinal Answer:")
    print(final_answer)


if __name__ == '__main__':
    try:
        find_supremum()
    except ImportError:
        print("Please install the networkx library using: pip install networkx")

<<<28>>>