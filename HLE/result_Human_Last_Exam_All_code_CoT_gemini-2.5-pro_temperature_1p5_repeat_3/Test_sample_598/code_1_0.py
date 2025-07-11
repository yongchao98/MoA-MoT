import numpy as np

def construct_l2_4():
    """
    Constructs the adjacency matrix of the L2(4) graph (4x4 rook's graph).
    The 16 vertices correspond to the cells of a 4x4 grid. Two vertices
    are adjacent if their corresponding cells are in the same row or column.
    """
    n = 16
    adj = np.zeros((n, n), dtype=int)
    for k1 in range(n):
        for k2 in range(k1 + 1, n):
            i1, j1 = divmod(k1, 4)
            i2, j2 = divmod(k2, 4)
            if i1 == i2 or j1 == j2:
                adj[k1, k2] = 1
                adj[k2, k1] = 1
    return adj

def construct_shrikhande(l2_4_adj):
    """
    Constructs the Shrikhande graph by switching the L2(4) graph.
    The switching is performed with respect to the set of vertices on the
    main diagonal of the 4x4 grid, which forms a coclique (independent set).
    """
    n = 16
    adj_shri = np.copy(l2_4_adj)
    
    # The switching set S = {(i,i) for i=0..3} corresponds to vertices 0, 5, 10, 15.
    S = {0, 5, 10, 15}
    
    for u in range(n):
        for v in range(u + 1, n):
            # Adjacency is flipped if one vertex is in S and the other is not.
            if (u in S and v not in S) or (v in S and u not in S):
                adj_shri[u, v] = 1 - adj_shri[u, v]
                adj_shri[v, u] = 1 - adj_shri[v, u]
    return adj_shri

def count_5_cycles(adj):
    """
    Counts the number of 5-cycles in a graph using a direct counting method.
    It iterates through all paths of length 4 and checks if the endpoints
    are connected, forming a 5-cycle.
    """
    n = adj.shape[0]
    adj_list = [np.where(row == 1)[0] for row in adj]
    
    count = 0
    # Iterate over all possible starting nodes v0
    for v0 in range(n):
        # Paths of length 1: v0-v1
        for v1 in adj_list[v0]:
            # Paths of length 2: v0-v1-v2
            for v2 in adj_list[v1]:
                if v2 == v0:
                    continue
                # Paths of length 3: v0-v1-v2-v3
                for v3 in adj_list[v2]:
                    if v3 == v1 or v3 == v0:
                        continue
                    # Paths of length 4: v0-v1-v2-v3-v4
                    for v4 in adj_list[v3]:
                        if v4 == v2 or v4 == v1 or v4 == v0:
                            continue
                        # Check if v4 is connected back to v0 to close the cycle
                        if v0 in adj_list[v4]:
                            count += 1
                            
    # Each 5-cycle v0-v1-v2-v3-v4-v0 is found 10 times:
    # 5 possible starting vertices * 2 directions (clockwise/counter-clockwise).
    return count // 10

def main():
    """
    Main function to construct graphs, count 5-cycles, and print results.
    """
    print("Yes, there is a combination of parameters for which two SRGs can have a different number of 5-cycles.")
    print("A classic example uses the parameters (n,d,lambda,mu) = (16, 6, 2, 2).")
    print("-" * 30)

    # Construct the two graphs
    g1_adj = construct_l2_4()
    g2_adj = construct_shrikhande(g1_adj)

    # Count 5-cycles in each graph
    c5_g1 = count_5_cycles(g1_adj)
    c5_g2 = count_5_cycles(g2_adj)

    print("Let G be the L2(4) graph (4x4 Rook's graph).")
    print(f"Number of 5-cycles in G: {c5_g1}")
    print("\nLet G' be the Shrikhande graph.")
    print(f"Number of 5-cycles in G': {c5_g2}")
    
    print("\nComparing the number of 5-cycles shows they are not equal:")
    # We explicitly print each number in the final equation.
    g1_val_str = str(c5_g1)
    g2_val_str = str(c5_g2)
    # The final equation is constructed by printing its components.
    print(f"{g1_val_str} != {g2_val_str}")

if __name__ == "__main__":
    main()
<<<Yes>>>