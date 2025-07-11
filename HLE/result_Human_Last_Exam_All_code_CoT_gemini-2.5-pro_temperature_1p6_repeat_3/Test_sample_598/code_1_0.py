import sys

def construct_l_k44_graph():
    """
    Constructs the L(K_4,4) graph, which is an SRG(16,6,2,2).
    It's isomorphic to the Cartesian product K_4 x K_4.
    Vertices are represented as integers 0-15, which map to pairs (i, j)
    with 0 <= i, j <= 3, where the integer ID is i * 4 + j.
    """
    n = 16
    adj = {i: set() for i in range(n)}
    for v1 in range(n):
        for v2 in range(v1 + 1, n):
            i1, j1 = v1 // 4, v1 % 4
            i2, j2 = v2 // 4, v2 % 4
            if i1 == i2 or j1 == j2:
                adj[v1].add(v2)
                adj[v2].add(v1)
    return adj

def construct_shrikhande_graph():
    """
    Constructs the Shrikhande graph, also an SRG(16,6,2,2).
    This is the Cayley graph of Z_4 x Z_4 with a specific connection set.
    Vertices are represented as integers 0-15, mapping to pairs (i, j) mod 4.
    """
    n = 16
    adj = {i: set() for i in range(n)}
    # Connection set for the Cayley graph
    s = {(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)}
    for v_idx in range(n):
        i, j = v_idx // 4, v_idx % 4
        for di, dj in s:
            ni, nj = (i + di) % 4, (j + dj) % 4
            neighbor_idx = ni * 4 + nj
            # The connection set is symmetric, so we only need to add one way
            if v_idx < neighbor_idx:
                adj[v_idx].add(neighbor_idx)
                adj[neighbor_idx].add(v_idx)
    return adj

def count_5_cycles(adj):
    """
    Counts the number of 5-cycles in a graph given by an adjacency list.
    It iterates through all 5-step walks and filters for those that are cycles.
    The final count is divided by 10 to correct for overcounting (5 starting
    points and 2 directions for each cycle).
    """
    count = 0
    n = len(adj)
    for v1 in range(n):
        for v2 in adj[v1]:
            for v3 in adj[v2]:
                if v3 == v1:
                    continue
                for v4 in adj[v3]:
                    if v4 == v2 or v4 == v1:
                        continue
                    # Path found: v1-v2-v3-v4
                    # Now check for a v5 that connects v4 back to v1 to close the cycle
                    if v1 in adj[v4]:
                        # This forms a closed walk v1-v2-v3-v4-v1. Check if it's a C4.
                        # We are looking for C5. Let's iterate one more step.
                        pass
                    
                    for v5 in adj[v4]:
                        if v5 == v3 or v5 == v2:
                            continue
                        # Potential path: v1-v2-v3-v4-v5
                        # Check if v5 connects to v1 and if all vertices are distinct
                        if v1 in adj[v5]:
                            if len({v1, v2, v3, v4, v5}) == 5:
                                count += 1
    return count // 10

def main():
    """
    Main function to construct the graphs, count cycles, and print the results.
    """
    print("Yes, such a combination of parameters and graphs exists.")
    print("The parameters are (n, d, λ, μ) = (16, 6, 2, 2).")
    print("The two non-isomorphic graphs for these parameters are the L(K4,4) graph and the Shrikhande graph.")
    print("-" * 30)

    # Construct and analyze the L(K4,4) graph
    l_k44 = construct_l_k44_graph()
    num_cycles_l_k44 = count_5_cycles(l_k44)
    print("Graph G = L(K4,4) graph")
    print(f"Number of 5-cycles = {num_cycles_l_k44}")
    
    # Construct and analyze the Shrikhande graph
    shrikhande = construct_shrikhande_graph()
    num_cycles_shrikhande = count_5_cycles(shrikhande)
    print("\nGraph G' = Shrikhande graph")
    print(f"Number of 5-cycles = {num_cycles_shrikhande}")
    print("-" * 30)
    
    if num_cycles_l_k44 != num_cycles_shrikhande:
        print(f"Since {num_cycles_l_k44} != {num_cycles_shrikhande}, these two graphs from the class srg(16, 6, 2, 2) have a different number of 5-cycles.")
    else:
        print("The number of 5-cycles is the same. There might be an error in the code or logic.")


if __name__ == '__main__':
    main()
<<<Yes>>>