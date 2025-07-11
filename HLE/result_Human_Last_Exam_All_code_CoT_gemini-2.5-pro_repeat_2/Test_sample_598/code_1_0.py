import networkx as nx

def count_5_cycles(G):
    """
    Counts the number of simple 5-cycles in a graph G.
    
    The method iterates through all paths of 4 vertices (v1-v2-v3-v4-v5)
    and checks if the last vertex (v5) is connected to the first (v1),
    while ensuring all vertices in the path are distinct.
    
    The final count is divided by 10 to correct for overcounting, as each
    cycle is found 5 times (once for each starting vertex) and 2 times
    (for each of the two directions of traversal).
    """
    count = 0
    # We convert nodes to a list to have a fixed order for iteration.
    nodes = list(G.nodes())
    for v1 in nodes:
        # Path: v1 -> v2
        for v2 in G.neighbors(v1):
            # Path: v1 -> v2 -> v3
            for v3 in G.neighbors(v2):
                if v3 == v1:
                    continue
                # Path: v1 -> v2 -> v3 -> v4
                for v4 in G.neighbors(v3):
                    if v4 == v1 or v4 == v2:
                        continue
                    # Path: v1 -> v2 -> v3 -> v4 -> v5
                    # We only need to check neighbors of v4 that are also neighbors of v1
                    # This is more efficient than a 5th nested loop.
                    common_neighbors = set(G.neighbors(v4)) & set(G.neighbors(v1))
                    for v5 in common_neighbors:
                        if v5 != v2 and v5 != v3:
                            count += 1
                            
    # Each 5-cycle v1-v2-v3-v4-v5-v1 is found when starting from each vertex.
    # For instance, v2-v3-v4-v5-v1-v2 is another found path.
    # The algorithm as written counts ordered tuples (v1,v2,v3,v4,v5)
    # where v1-v2-...-v5-v1 is a cycle. There are 2*5=10 such tuples for each cycle.
    # However, my optimized loop counts pairs (v1, cycle) which gives 5 counts per cycle.
    # And the inner part `for v5 in common_neighbors` counts (v1,v2,v3,v4,v5) and (v1,v5,v4,v3,v2)
    # so we divide by 2*5 = 10. Let's trace it.
    # For a cycle (a,b,c,d,e), my code will count:
    # (a,b,c,d,e), (a,e,d,c,b), (b,c,d,e,a), (b,a,e,d,c), etc.
    # Total count will be 2*5 times the number of cycles.
    return count // 10

# 1. Define the parameters
n = 16
d = 6
lambda_param = 2
mu_param = 2

print(f"The chosen SRG class is srg(n,d,lambda,mu) = srg({n}, {d}, {lambda_param}, {mu_param}).")
print("-" * 60)

# 2. Generate the two graphs
# Graph 1: The L2(4) graph, or 4x4 Rook's graph
G_rook = nx.rook_graph(4)

# Graph 2: The Shrikhande graph
G_shrikhande = nx.shrikhande_graph()

# 3. Count 5-cycles in each graph
cycles_rook = count_5_cycles(G_rook)
cycles_shrikhande = count_5_cycles(G_shrikhande)

# 4. Print the results
print("Graph 1: 4x4 Rook's graph (L2(4) graph)")
print(f"  - n (vertices): {G_rook.number_of_nodes()}")
print(f"  - d (degree): {G_rook.degree(0)}")
print(f"  - Number of 5-cycles: {cycles_rook}")
print()
print("Graph 2: Shrikhande graph")
print(f"  - n (vertices): {G_shrikhande.number_of_nodes()}")
print(f"  - d (degree): {G_shrikhande.degree(0)}")
print(f"  - Number of 5-cycles: {cycles_shrikhande}")
print("-" * 60)

if cycles_rook != cycles_shrikhande:
    print("Conclusion: The two graphs are both in srg(16, 6, 2, 2) but have a different number of 5-cycles.")
else:
    print("Conclusion: The two graphs have the same number of 5-cycles.")
