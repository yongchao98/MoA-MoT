import networkx as nx

def create_rooks_graph():
    """
    Constructs the 4x4 Rook's graph, which is an srg(16, 6, 2, 2).
    Vertices are integers 0-15, representing cells of a 4x4 grid.
    v = 4*r + c.
    Two vertices are adjacent if they are in the same row or column.
    """
    G = nx.Graph()
    n_size = 4
    nodes = range(n_size * n_size)
    G.add_nodes_from(nodes)
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            r1, c1 = divmod(i, n_size)
            r2, c2 = divmod(j, n_size)
            if r1 == r2 or c1 == c2:
                G.add_edge(i, j)
    return G

def create_shrikhande_graph():
    """
    Constructs the Shrikhande graph, which is also an srg(16, 6, 2, 2).
    It is constructed as a Cayley graph on Z_4 x Z_4.
    Vertices are tuples (r, c) where r, c are in {0, 1, 2, 3}.
    """
    G = nx.Graph()
    n_size = 4
    nodes = [(r, c) for r in range(n_size) for c in range(n_size)]
    G.add_nodes_from(nodes)
    
    # Connection set S for the Shrikhande graph
    S = [(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)]
    
    for r1 in range(n_size):
        for c1 in range(n_size):
            for dr, dc in S:
                r2 = (r1 + dr) % n_size
                c2 = (c1 + dc) % n_size
                G.add_edge((r1, c1), (r2, c2))
    return G

def count_cycles(G, k):
    """
    Counts the number of simple cycles of length k in a graph G.
    It iterates through all simple paths of length k-1 and checks
    if the endpoints are connected.
    The total count is divided by 2*k to correct for double counting
    (each direction and each starting point).
    """
    count = 0
    for start_node in G.nodes():
        # Stack stores tuples of (current_node, path_list) for DFS
        stack = [(neighbor, [start_node, neighbor]) for neighbor in G.neighbors(start_node)]
        
        while stack:
            current_node, path = stack.pop()
            
            if len(path) == k:
                # Path of length k-1 found, check for closing edge
                if G.has_edge(path[-1], start_node):
                    count += 1
                continue
            
            # Extend the path with new neighbors
            for neighbor in G.neighbors(current_node):
                if neighbor not in path:
                    new_path = path + [neighbor]
                    stack.append((neighbor, new_path))
                    
    return count // (2 * k)

def main():
    """
    Main function to perform the graph creation, cycle counting, and result printing.
    """
    n, d, l, m = 16, 6, 2, 2
    k_cycle = 5

    print(f"Investigating strongly regular graphs with parameters (n, d, λ, μ) = ({n}, {d}, {l}, {m}).")
    print("-" * 60)

    # 4x4 Rook's Graph
    g1 = create_rooks_graph()
    c5_g1 = count_cycles(g1, k_cycle)
    print(f"Graph 1: The 4x4 Rook's Graph (L2(4))")
    print(f"Number of {k_cycle}-cycles: {c5_g1}")
    print("-" * 60)

    # Shrikhande Graph
    g2 = create_shrikhande_graph()
    c5_g2 = count_cycles(g2, k_cycle)
    print(f"Graph 2: The Shrikhande Graph")
    print(f"Number of {k_cycle}-cycles: {c5_g2}")
    print("-" * 60)

    if c5_g1 != c5_g2:
        print("Yes, there exists a combination of parameters (n,d,λ,μ) and a pair of graphs G, G'")
        print("in srg(n,d,λ,μ) such that G and G' do not have the same number of 5-cycles.")
        print(f"\nFor parameters ({n}, {d}, {l}, {m}):")
        print(f"The 4x4 Rook's graph has {c5_g1} {k_cycle}-cycles.")
        print(f"The Shrikhande graph has {c5_g2} {k_cycle}-cycles.")
    else:
        print("For the tested pair of graphs, the number of 5-cycles is the same.")

if __name__ == '__main__':
    main()
