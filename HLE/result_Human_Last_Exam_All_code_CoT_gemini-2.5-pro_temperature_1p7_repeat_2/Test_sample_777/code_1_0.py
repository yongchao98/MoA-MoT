import networkx as nx
from itertools import combinations

def solve_disjoint_cycles(graph_adj, k):
    """
    Solves the DisjointCycles problem using a backtracking approach.
    This is not an FPT algorithm and will be very slow for large graphs.

    Args:
        graph_adj (dict): The graph represented as an adjacency list.
        k (int): The number of disjoint cycles to find, and the minimum length for each.

    Returns:
        A tuple (result, cycles) where result is 1 or 0, and cycles is the list
        of found cycles if result is 1.
    """
    if k <= 0:
        return 1, []

    G = nx.Graph(graph_adj)
    
    # Step 1: Find all simple cycles and filter by length
    long_cycles = []
    # networkx.simple_cycles finds all elementary cycles in a directed graph
    # For an undirected graph, we can convert to directed considering each edge in both directions
    # A simple cycle of length L in G corresponds to a simple cycle of length L in D.
    D = G.to_directed()
    for cycle in nx.simple_cycles(D):
        # simple_cycles in networkx might produce duplicates for undirected graphs (e.g. A->B->C->A and A->C->B->A)
        # We normalize by using a frozenset.
        if len(cycle) >= k:
            long_cycles.append(frozenset(cycle))
    
    # Remove duplicate cycles
    unique_long_cycles = list(set(long_cycles))

    # Step 2: Search for k disjoint cycles among the candidates
    # We can frame this as finding a k-clique in the disjointness graph, but a simple
    # combination check is easier to write for a demonstration.
    if len(unique_long_cycles) < k:
        print(f"Found only {len(unique_long_cycles)} cycles of length >= {k}. Cannot find {k} disjoint cycles.")
        return 0, []

    for combo in combinations(unique_long_cycles, k):
        is_disjoint_set = True
        # Check for pairwise disjointness in the combination
        for i in range(k):
            for j in range(i + 1, k):
                if not combo[i].isdisjoint(combo[j]):
                    is_disjoint_set = False
                    break
            if not is_disjoint_set:
                break
        
        if is_disjoint_set:
            # Found a valid set
            solution = [list(c) for c in combo]
            print(f"Success! Found {k} vertex-disjoint cycles of length >= {k}:")
            for i, cycle in enumerate(solution):
                print(f"  Cycle {i+1} (length {len(cycle)}): {cycle}")
            return 1, solution
            
    print(f"Failure. Could not find {k} vertex-disjoint cycles of length >= {k}.")
    return 0, []


# --- Example Usage ---

# A graph with a 5-cycle {0,1,2,3,4} and a 4-cycle {5,6,7,8}
graph = {
    0: [1, 4], 1: [0, 2], 2: [1, 3], 3: [2, 4], 4: [0, 3],
    5: [6, 8], 6: [5, 7], 7: [6, 8], 8: [5, 7]
}

# Case 1: k=2. We need >= 2 cycles of length >= 2. The graph has them.
print("--- Running for k=2 ---")
k_val = 2
result_code, _ = solve_disjoint_cycles(graph, k_val)
print(f"Problem: DisjointCycles(G, k={k_val})")
print(f"Output: {result_code}")
print("-" * 20)


# Case 2: k=3. We need >= 3 cycles of length >= 3. The graph only has two.
print("--- Running for k=3 ---")
k_val = 3
result_code, _ = solve_disjoint_cycles(graph, k_val)
print(f"Problem: DisjointCycles(G, k={k_val})")
print(f"Output: {result_code}")
print("-" * 20)

# Case 3: k=5. We need >= 5 cycles of length >= 5. The graph only has one such cycle.
print("--- Running for k=5 ---")
k_val = 5
result_code, _ = solve_disjoint_cycles(graph, k_val)
print(f"Problem: DisjointCycles(G, k={k_val})")
print(f"Output: {result_code}")
print("-" * 20)
