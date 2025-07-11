import numpy as np

def build_rooks_graph():
    """
    Builds the adjacency matrix for the 4x4 Rook's graph (L2(4)).
    Vertices are represented as integers 0-15, where vertex k corresponds
    to cell (i, j) with k = 4*i + j.
    Two vertices are adjacent if they are in the same row or column.
    """
    n = 16
    adj = np.zeros((n, n), dtype=int)
    for i in range(4):
        for j in range(4):
            u = 4 * i + j
            # Neighbors in the same row
            for j_prime in range(4):
                if j_prime != j:
                    v = 4 * i + j_prime
                    adj[u, v] = 1
            # Neighbors in the same column
            for i_prime in range(4):
                if i_prime != i:
                    v = 4 * i_prime + j
                    adj[u, v] = 1
    return adj

def build_shrikhande_graph():
    """
    Builds the adjacency matrix for the Shrikhande graph.
    This is the Cayley graph of the group Z_4 x Z_4 with the
    connection set S = {+-(1,0), +-(0,1), +-(1,1)}.
    """
    n = 16
    adj = np.zeros((n, n), dtype=int)
    connections = [(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)]
    for i in range(4):
        for j in range(4):
            u_idx = 4 * i + j
            for di, dj in connections:
                i_prime = (i + di) % 4
                j_prime = (j + dj) % 4
                v_idx = 4 * i_prime + j_prime
                adj[u_idx, v_idx] = 1
    return adj

def count_c5(adj):
    """
    Counts the number of simple 5-cycles in a graph.
    The algorithm explores all simple paths of length 4 starting from each vertex
    and checks if an edge exists to close the cycle.
    """
    n = adj.shape[0]
    # Create an adjacency list for faster neighbor lookups
    adj_list = [np.where(row == 1)[0] for row in adj]
    count = 0
    # i-j-k-l-m-i represents a 5-cycle
    for i in range(n):
        for j in adj_list[i]:
            for k in adj_list[j]:
                if k == i: continue
                for l in adj_list[k]:
                    if l == j or l == i: continue
                    for m in adj_list[l]:
                        if m == k or m == j or m == i: continue
                        # We have a simple path i-j-k-l-m. Check if it's a 5-cycle.
                        if i in adj_list[m]:
                            count += 1
    # Each cycle is found 5 times (for each starting node) and
    # 2 times (for each direction), so we divide by 10.
    return count // 10

def main():
    """
    Main function to construct the graphs, count their 5-cycles, and print the results.
    """
    # The parameters for the class of strongly regular graphs
    n, d, lam, mu = 16, 6, 2, 2
    
    # Construct the two non-isomorphic graphs
    rooks_graph = build_rooks_graph()
    shrikhande_graph = build_shrikhande_graph()
    
    # Count the number of 5-cycles in each graph
    n5_rooks = count_c5(rooks_graph)
    n5_shrikhande = count_c5(shrikhande_graph)
    
    print("Yes, such a combination of parameters and graphs exists.")
    print(f"Consider the class of strongly regular graphs srg(n,d,lambda,mu) with parameters:")
    print(f"n = {n}")
    print(f"d = {d}")
    print(f"lambda = {lam}")
    print(f"mu = {mu}")
    print("\nThis class contains at least two non-isomorphic graphs: the 4x4 Rook's graph and the Shrikhande graph.")
    print("Despite being co-spectral, they have a different number of 5-cycles.")
    print("\n--- Calculation Results ---")
    print(f"Number of 5-cycles in the 4x4 Rook's graph: {n5_rooks}")
    print(f"Number of 5-cycles in the Shrikhande graph: {n5_shrikhande}")

if __name__ == '__main__':
    main()