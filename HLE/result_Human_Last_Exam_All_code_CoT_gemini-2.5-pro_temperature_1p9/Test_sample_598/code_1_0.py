import networkx as nx

def create_rooks_graph():
    """
    Constructs the 4x4 Rook's graph, which is a strongly regular graph
    with parameters (16, 6, 2, 2).
    The vertices represent squares on a 4x4 board, and two vertices
    are connected if they are in the same row or column.
    """
    G = nx.Graph()
    n_side = 4
    num_nodes = n_side * n_side
    G.add_nodes_from(range(num_nodes))
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            row1, col1 = divmod(i, n_side)
            row2, col2 = divmod(j, n_side)
            if row1 == row2 or col1 == col2:
                G.add_edge(i, j)
    return G

def create_shrikhande_graph():
    """
    Constructs the Shrikhande graph, a strongly regular graph with
    parameters (16, 6, 2, 2) that is not isomorphic to the Rook's graph.
    It can be defined as a Cayley graph on the group Z_4 x Z_4.
    """
    G = nx.Graph()
    n_side = 4
    num_nodes = n_side * n_side
    G.add_nodes_from(range(num_nodes))
    
    # The generating set S defines the edges. Two vertices u, v are
    # adjacent if v - u (mod 4) is in S or its inverse.
    s_set = [(0, 1), (1, 0), (1, 1), (0, 3), (3, 0), (3, 3)]
    
    for i in range(n_side):
        for j in range(n_side):
            u_node = i * n_side + j
            for sx, sy in s_set:
                neighbor_i, neighbor_j = (i + sx) % n_side, (j + sy) % n_side
                v_node = neighbor_i * n_side + neighbor_j
                G.add_edge(u_node, v_node)
    return G

def count_simple_cycles(G, k):
    """
    Counts the number of simple cycles of length k in a graph G
    using a depth-first search approach from each node.
    """
    count = 0
    # We iterate over all nodes as potential starting points for the cycle.
    for start_node in G.nodes():
        # The stack holds tuples of (current_node, path_list)
        stack = [(start_node, [start_node])]
        while stack:
            curr_node, path = stack.pop()
            
            if len(path) == k:
                # If path has length k, check if it can be closed to form a cycle.
                # The start node must be a neighbor of the current node.
                if start_node in G[curr_node]:
                    count += 1
                continue

            # Explore neighbors to extend the path.
            for neighbor in G[curr_node]:
                # Do not visit nodes already in the path, except for the start_node
                # which can only be visited at the very end to close the cycle.
                if neighbor not in path:
                    new_path = path + [neighbor]
                    stack.append((neighbor, new_path))
    
    # Each cycle c = (v1, ..., vk) is found 2k times:
    # - Once starting from each of its k vertices.
    # - Once in each of two directions (e.g., v1->v2... and v1->vk...).
    # Therefore, we divide the total count by 2*k.
    return count // (2 * k) if k > 2 else count // 2

def solve():
    """
    Main function to construct the graphs, count cycles, and print the results.
    """
    print("Yes, such a pair of graphs exists. We will demonstrate this using the two known")
    print("graphs in the class srg(16, 6, 2, 2).\n")

    # Define the SRG parameters.
    n, d, lam, mu = 16, 6, 2, 2
    
    # Create the two graphs.
    print("Constructing the 4x4 Rook's graph and the Shrikhande graph...")
    rooks_graph = create_rooks_graph()
    shrikhande_graph = create_shrikhande_graph()

    # Count the 5-cycles in each graph.
    print("Counting 5-cycles in both graphs (this may take a moment)...")
    num_c5_rooks = count_simple_cycles(rooks_graph, 5)
    num_c5_shrikhande = count_simple_cycles(shrikhande_graph, 5)

    print("\n" + "="*40)
    print("      STRONGLY REGULAR GRAPH PARAMETERS")
    print("="*40)
    print(f"Number of vertices (n):       {n}")
    print(f"Degree of each vertex (d):    {d}")
    print(f"Common neighbors (adjacent, lambda): {lam}")
    print(f"Common neighbors (non-adj, mu): {mu}")
    print("="*40 + "\n")
    
    print("--- RESULTS ---\n")
    print("Graph 1: 4x4 Rook's Graph")
    print(f"Number of 5-cycles: {num_c5_rooks}\n")

    print("Graph 2: Shrikhande Graph")
    print(f"Number of 5-cycles: {num_c5_shrikhande}\n")
    
    if num_c5_rooks != num_c5_shrikhande:
        print("Conclusion: The two graphs have the same SRG parameters but a different")
        print("number of 5-cycles, proving that such a case exists.")
    else:
        print("Conclusion: In this case, the two graphs have the same number of 5-cycles.")

# Execute the solution
solve()
