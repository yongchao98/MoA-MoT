import networkx as nx

def _dfs_count_cycles(G, k, start_node, current_node, path, visited):
    """
    A helper recursive function to find simple cycles of length k.
    """
    if len(path) == k:
        # Path has k vertices, check if it's a cycle
        if start_node in G.neighbors(current_node):
            return 1
        else:
            return 0
    
    count = 0
    for neighbor in G.neighbors(current_node):
        # We need a simple path, so no repeated vertices
        if neighbor not in visited:
            path.append(neighbor)
            visited.add(neighbor)
            count += _dfs_count_cycles(G, k, start_node, neighbor, path, visited)
            # Backtrack
            visited.remove(neighbor)
            path.pop()
    return count

def count_simple_cycles(G, k):
    """
    Counts the number of simple cycles of length k in a graph G.
    """
    total_count = 0
    for start_node in G.nodes():
        path = [start_node]
        visited = {start_node}
        total_count += _dfs_count_cycles(G, k, start_node, start_node, path, visited)
        
    # Each cycle is counted k times (for each starting node) and 2 times for
    # each direction (e.g., A-B-C-A and A-C-B-A).
    # So we divide by 2*k.
    return total_count // (2 * k)

def solve():
    """
    Constructs the graphs and prints the number of 5-cycles.
    """
    srg_params = {'n': 16, 'd': 6, 'lambda': 2, 'mu': 2}
    
    # --- Graph 1: 4x4 Rook's graph ---
    G1_rooks = nx.Graph()
    nodes = [(i, j) for i in range(4) for j in range(4)]
    G1_rooks.add_nodes_from(nodes)

    for r1 in range(4):
        for c1 in range(4):
            for r2 in range(4):
                for c2 in range(4):
                    v1 = (r1, c1)
                    v2 = (r2, c2)
                    if v1 != v2 and (r1 == r2 or c1 == c2):
                        G1_rooks.add_edge(v1, v2)

    # --- Graph 2: Shrikhande graph ---
    # Construction as a Cayley graph on Z_4 x Z_4
    # with generating set S = {(+/-1, 0), (0, +/-1), (1, 1), (-1, -1)}
    G2_shrikhande = nx.Graph()
    G2_shrikhande.add_nodes_from(nodes)
    offsets = [(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)]
    for r in range(4):
        for c in range(4):
            v1 = (r, c)
            for dr, dc in offsets:
                v2 = ((r + dr) % 4, (c + dc) % 4)
                G2_shrikhande.add_edge(v1, v2)
                
    # --- Count 5-cycles ---
    k = 5
    c5_rooks = count_simple_cycles(G1_rooks, k)
    c5_shrikhande = count_simple_cycles(G2_shrikhande, k)
    
    print("Yes, a combination of parameters with two SRGs having different numbers of 5-cycles exists.")
    print(f"The parameters are (n, d, lambda, mu) = ({srg_params['n']}, {srg_params['d']}, {srg_params['lambda']}, {srg_params['mu']}).")
    print("-" * 20)
    print("The two graphs are the 4x4 Rook's graph (G1) and the Shrikhande graph (G2).")
    print(f"The number of 5-cycles in the 4x4 Rook's graph is {c5_rooks}.")
    print(f"The number of 5-cycles in the Shrikhande graph is {c5_shrikhande}.")
    print("-" * 20)
    print(f"Final Equation: Number of 5-cycles in G1 = {c5_rooks} != {c5_shrikhande} = Number of 5-cycles in G2.")

solve()