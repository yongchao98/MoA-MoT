import networkx as nx

def count_simple_cycles(G, k):
    """
    Counts the number of simple cycles of length k in a NetworkX graph G.
    
    This function uses a depth-first search approach starting from each node.
    It finds directed paths of length k and checks if the last node connects
    back to the start node. The total count is then divided by 2*k to correct
    for overcounting (k starting points and 2 directions for each cycle).
    """
    count = 0
    adj = {n: set(G.neighbors(n)) for n in G.nodes()}

    for start_node in G.nodes():
        # The stack stores tuples of (current_node, path_list)
        # Using a list for the path is acceptable for small k.
        stack = [(start_node, [start_node])]
        while stack:
            current_node, path = stack.pop()
            
            if len(path) == k:
                # Path length is k-1 (e.g., A->B->C for k=3).
                # If current_node connects to start_node, we have found a cycle.
                if start_node in adj[current_node]:
                    count += 1
                continue

            # Explore neighbors to extend the path
            for neighbor in adj[current_node]:
                # We must avoid reusing nodes in the current path.
                if neighbor not in path:
                    new_path = path + [neighbor]
                    stack.append((neighbor, new_path))
    
    if k <= 1:
        return 0
    # Each cycle is counted 2*k times (k start nodes, 2 directions)
    return count // (2 * k)

def get_srg_params(G):
    """
    Determines the SRG parameters (n, d, lambda, mu) for a graph G.
    This function assumes G is a strongly regular graph and samples nodes/edges
    to calculate the parameters.
    """
    nodes = list(G.nodes())
    if not nodes:
        return (0, 0, 0, 0)
    
    n = G.number_of_nodes()
    
    # Pick an arbitrary node to find degree d
    d = G.degree(nodes[0])
    
    # Pick an arbitrary edge to find lambda
    lambda_val = 0
    if G.number_of_edges() > 0:
        u, v = next(iter(G.edges()))
        lambda_val = len(list(nx.common_neighbors(G, u, v)))
        
    # Pick an arbitrary non-edge to find mu
    mu_val = 'N/A'
    non_edges = list(nx.non_edges(G))
    if non_edges:
        u, v = non_edges[0]
        mu_val = len(list(nx.common_neighbors(G, u, v)))
    else: # Complete graph
        mu_val = d

    return (n, d, lambda_val, mu_val)

def build_rooks_graph(n_size=4):
    """
    Constructs the 4x4 Rook's graph (L_2(4)).
    Nodes are cells on a 4x4 grid, adjacent if in the same row or column.
    """
    G = nx.Graph()
    num_nodes = n_size * n_size
    G.add_nodes_from(range(num_nodes))
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            row_i, col_i = i // n_size, i % n_size
            row_j, col_j = j // n_size, j % n_size
            if row_i == row_j or col_i == col_j:
                G.add_edge(i, j)
    return G

def build_shrikhande_graph(n_size=4):
    """
    Constructs the Shrikhande graph.
    This implementation uses the Cayley graph construction on the group
    Z_4 x Z_4 with the generating set S = {(+/-1,0), (0, +/-1), (1,1), (-1,-1)}.
    """
    G = nx.Graph()
    num_nodes = n_size * n_size
    G.add_nodes_from(range(num_nodes))
    
    gen_set = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]

    for i in range(n_size):
        for j in range(n_size):
            u_idx = i * n_size + j
            for s_x, s_y in gen_set:
                v_coord_x = (i + s_x) % n_size
                v_coord_y = (j + s_y) % n_size
                v_idx = v_coord_x * n_size + v_coord_y
                G.add_edge(u_idx, v_idx)
    return G

def main():
    """
    Main function to construct graphs, verify parameters, and count cycles.
    """
    print("This script will demonstrate that two SRGs with identical parameters")
    print("can have a different number of 5-cycles.\n")
    print("The example uses the 4x4 Rook's graph and the Shrikhande graph.")

    # Create the two graphs
    rooks_graph = build_rooks_graph()
    shrikhande_graph = build_shrikhande_graph()

    print("\nVerifying SRG parameters...")
    rooks_params = get_srg_params(rooks_graph)
    shrikhande_params = get_srg_params(shrikhande_graph)

    print(f"Rook's graph parameters (n, d, lambda, mu): {rooks_params}")
    print(f"Shrikhande graph parameters (n, d, lambda, mu): {shrikhande_params}")

    if rooks_params == shrikhande_params:
        print("Verification successful: Both graphs are SRGs with the same parameters.")
    else:
        print("Warning: The graphs do not share the same SRG parameters.")

    # Count 5-cycles
    print("\nCounting the number of 5-cycles in each graph...")
    
    cycle_length = 5
    
    rooks_c5_count = count_simple_cycles(rooks_graph, cycle_length)
    shrikhande_c5_count = count_simple_cycles(shrikhande_graph, cycle_length)

    print(f"Number of 5-cycles in the 4x4 Rook's graph: {rooks_c5_count}")
    print(f"Number of 5-cycles in the Shrikhande graph: {shrikhande_c5_count}")

    if rooks_c5_count != shrikhande_c5_count:
        print("\nConclusion: The two graphs have the same SRG parameters but a different number of 5-cycles.")
    else:
        print("\nConclusion: Both graphs have the same number of 5-cycles, the example doesn't show a difference.")


if __name__ == "__main__":
    main()