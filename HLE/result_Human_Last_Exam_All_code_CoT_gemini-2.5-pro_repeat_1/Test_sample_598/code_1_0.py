import numpy as np

def create_rook_graph_4x4():
    """
    Constructs the adjacency matrix of the L(K_4,4) graph (4x4 Rook's graph).
    This is an srg(16, 6, 2, 2).
    Vertices are cells (i, j) of a 4x4 grid, 0 <= i,j <= 3.
    Two vertices are adjacent if they are in the same row or column.
    """
    n = 16
    adj = np.zeros((n, n), dtype=int)
    for v1 in range(n):
        for v2 in range(v1 + 1, n):
            i1, j1 = v1 // 4, v1 % 4
            i2, j2 = v2 // 4, v2 % 4
            if i1 == i2 or j1 == j2:
                adj[v1, v2] = adj[v2, v1] = 1
    return adj

def create_shrikhande_graph():
    """
    Constructs the adjacency matrix of the Shrikhande graph.
    This is also an srg(16, 6, 2, 2).
    It's the Cayley graph of Z_4 x Z_4 with connection set
    S = { (1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1) }.
    """
    n = 16
    adj = np.zeros((n, n), dtype=int)
    # In Z_4, -1 is 3.
    gen_set = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
    
    for v1 in range(n):
        i1, j1 = v1 // 4, v1 % 4
        for s1, s2 in gen_set:
            i2 = (i1 + s1) % 4
            j2 = (j1 + s2) % 4
            v2 = 4 * i2 + j2
            adj[v1, v2] = 1
    return adj

def count_5_cycles(adj):
    """
    Counts the number of 5-cycles (C5) in a graph given its adjacency matrix.
    
    The algorithm iterates through all paths of length 4 (v1-v2-v3-v4-v5)
    and checks if the last vertex v5 is connected to the first vertex v1,
    forming a cycle. It ensures all vertices in the path are distinct.
    
    Each cycle is found 10 times (5 starting points * 2 directions), so the
    final count is divided by 10.
    """
    n = len(adj)
    count = 0
    
    # Create neighbor lists for faster access
    neighbors = [np.where(row == 1)[0] for row in adj]

    for v1 in range(n):
        for v2 in neighbors[v1]:
            for v3 in neighbors[v2]:
                if v3 == v1:
                    continue
                for v4 in neighbors[v3]:
                    if v4 == v1 or v4 == v2:
                        continue
                    for v5 in neighbors[v4]:
                        if v5 == v1 or v5 == v2 or v5 == v3:
                            continue
                        # We have a distinct path v1-v2-v3-v4-v5.
                        # Check for the closing edge (v5, v1).
                        if v1 in neighbors[v5]:
                            count += 1
                            
    return count // 10

def main():
    """
    Main function to construct the graphs, count their 5-cycles, and print the results.
    """
    params = {'n': 16, 'd': 6, 'lambda': 2, 'mu': 2}
    
    print(f"Investigating strongly regular graphs with parameters (n,d,lambda,mu) = ({params['n']},{params['d']},{params['lambda']},{params['mu']}).")
    print("-" * 70)

    # Rook's graph
    adj_rook = create_rook_graph_4x4()
    c5_rook = count_5_cycles(adj_rook)
    print("Graph 1: The 4x4 Rook's graph (L(K_4,4))")
    print(f"  - Number of vertices (n): {params['n']}")
    print(f"  - Degree (d): {params['d']}")
    print(f"  - Common neighbors for adjacent vertices (lambda): {params['lambda']}")
    print(f"  - Common neighbors for non-adjacent vertices (mu): {params['mu']}")
    print(f"  - Number of 5-cycles: {c5_rook}")
    print("-" * 70)

    # Shrikhande graph
    adj_shrikhande = create_shrikhande_graph()
    c5_shrikhande = count_5_cycles(adj_shrikhande)
    print("Graph 2: The Shrikhande graph")
    print(f"  - Number of vertices (n): {params['n']}")
    print(f"  - Degree (d): {params['d']}")
    print(f"  - Common neighbors for adjacent vertices (lambda): {params['lambda']}")
    print(f"  - Common neighbors for non-adjacent vertices (mu): {params['mu']}")
    print(f"  - Number of 5-cycles: {c5_shrikhande}")
    print("-" * 70)
    
    if c5_rook != c5_shrikhande:
        print("Conclusion: Yes, there are SRGs with the same parameters but a different number of 5-cycles.")
    else:
        print("Conclusion: No, these two graphs have the same number of 5-cycles.")

if __name__ == '__main__':
    main()
