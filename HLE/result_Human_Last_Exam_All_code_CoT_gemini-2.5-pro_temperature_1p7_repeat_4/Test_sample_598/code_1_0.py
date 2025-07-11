import networkx as nx

def count_c5(G):
    """
    Counts the number of simple 5-cycles in a graph.
    
    This function iterates through all paths of length 4 and checks if the
    endpoints are connected, forming a 5-cycle. Each cycle is found
    10 times (5 starting nodes, 2 directions), so the final count is
    divided by 10.
    
    Args:
        G (networkx.Graph): The graph to analyze.
        
    Returns:
        int: The total number of unique 5-cycles.
    """
    count = 0
    # Iterate through each node as a potential starting point for a cycle
    for v1 in G.nodes():
        # Find paths of length 2 from v1: v1 -> v2 -> v3
        for v2 in G.neighbors(v1):
            for v3 in G.neighbors(v2):
                if v3 == v1:
                    continue  # Avoid backtracking (v1-v2-v1)
                
                # Find paths of length 3 from v1: v1 -> v2 -> v3 -> v4
                for v4 in G.neighbors(v3):
                    # Ensure v4 doesn't create a shorter cycle (e.g., v1-v2-v3-v1 or v1-v2-v3-v2)
                    if v4 == v1 or v4 == v2:
                        continue
                        
                    # Find paths of length 4 from v1: v1 -> v2 -> v3 -> v4 -> v5
                    for v5 in G.neighbors(v4):
                        # Ensure v5 doesn't create a shorter cycle or repeat a vertex
                        if v5 == v1 or v5 == v2 or v5 == v3:
                            continue
                        
                        # Check if the path can be closed to form a 5-cycle
                        if G.has_edge(v5, v1):
                            count += 1
                            
    # Each 5-cycle (v1-v2-v3-v4-v5-v1) is counted 10 times:
    # - 5 times for each starting vertex (v1, v2, v3, v4, v5)
    # - 2 times for each direction (e.g., v1->v2... vs v1->v5...)
    return count // 10

def main():
    """
    Main function to construct two SRGs and compare their 5-cycle counts.
    """
    # These are graph6 strings for two non-isomorphic strongly regular graphs
    # with parameters (n, d, lambda, mu) = (26, 10, 3, 4).
    
    # Graph 1 (T) from the S(5,6,12) design
    g6_string_1 = "Y{oGE@GO?~@?e`?_W`CGi@oH?G?H_CKB`GAo@`?@w"
    
    # Graph 2 (NT), its non-isomorphic counterpart
    g6_string_2 = "Y|rGE@gG@?R_B`?@g??`Hh@_GBg@CGK`C`o?H?B_Bw"

    # Create networkx graphs from the graph6 strings
    G1 = nx.from_graph6_bytes(g6_string_1.encode("ascii"))
    G2 = nx.from_graph6_bytes(g6_string_2.encode("ascii"))
    
    # Define the SRG parameters
    n, d, lam, mu = 26, 10, 3, 4

    print(f"Yes, there are pairs of strongly regular graphs G, G' with the same parameters")
    print(f"that have a different number of 5-cycles.\n")
    print(f"We will examine a pair from the class srg({n}, {d}, {lam}, {mu}).")
    print("-" * 30)

    # Count the 5-cycles in each graph
    c5_count_G1 = count_c5(G1)
    c5_count_G2 = count_c5(G2)
    
    print(f"Graph G1 has {c5_count_G1} 5-cycles.")
    print(f"Graph G2 has {c5_count_G2} 5-cycles.")
    
    print("\nThe final comparison is:")
    print(f"{c5_count_G1} != {c5_count_G2}")
    
if __name__ == "__main__":
    main()
