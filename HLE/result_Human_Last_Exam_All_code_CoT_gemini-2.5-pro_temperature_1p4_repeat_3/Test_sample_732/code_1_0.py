import collections

def solve_max_c4_free_edges():
    """
    This function calculates the maximum number of edges in a C4-free graph on 8 vertices.
    It combines a theoretical upper bound analysis with a computational search for the exact number.
    """
    n = 8
    
    print("Step 1: Understanding the problem and deriving a theoretical upper bound.")
    print("The problem asks for the maximum number of edges in a simple graph with 8 vertices such that it contains no quadrilaterals (C4).")
    print("A C4-free graph has the property that any two distinct vertices share at most one common neighbor.")
    print("This property leads to an inequality involving the degrees of the vertices ('d_i').")
    print("The inequality is: sum(C(d_i, 2) for i in 1..n) <= C(n, 2)")
    
    c_n_2 = n * (n - 1) // 2
    print(f"For n = {n}, C(n, 2) = C({n}, 2) = {c_n_2}.")
    print(f"So, the sum of C(d_i, 2) for all 8 vertices must be less than or equal to {c_n_2}.")
    
    print("\nThis inequality gives an upper bound for the number of edges, 'm'.")
    print("Based on this, the number of edges 'm' can be at most 12.")
    print("A graph with 8 vertices and 12 edges must be 3-regular (all vertices have degree 3).")
    print("It is known that all 3-regular graphs on 8 vertices contain at least one C4.")
    print("Therefore, the maximum number of edges must be less than 12.")
    
    print("\nStep 2: Performing a computational search for the exact maximum.")
    print("We use a backtracking algorithm to build a C4-free graph, adding edges one by one to find the largest possible graph.")

    # --- Backtracking search implementation ---
    # Use a list for max_edges_found to be mutable inside the nested function
    max_edges_found = [0] 
    
    # The search is optimized: a new C4 must involve the newly added edge (u, v).
    # This only happens if u and v already share a common neighbor before the edge is added.
    def fast_search(potential_edges, adj, start_index):
        # Update the maximum if a larger graph is found
        num_current_edges = sum(len(neighbors) for neighbors in adj.values()) // 2
        if num_current_edges > max_edges_found[0]:
            max_edges_found[0] = num_current_edges

        # Pruning: if we can't possibly beat the current max, stop.
        if num_current_edges + (len(potential_edges) - start_index) <= max_edges_found[0]:
             return

        for i in range(start_index, len(potential_edges)):
            u, v = potential_edges[i]
            
            # Check if adding edge (u, v) creates a C4
            if not adj[u].intersection(adj[v]):
                # Safe to add the edge
                adj[u].add(v)
                adj[v].add(u)
                
                fast_search(potential_edges, adj, i + 1)
                
                # Backtrack: remove the edge
                adj[u].remove(v)
                adj[v].remove(u)

    # Generate all possible edges for a simple graph on n vertices
    all_possible_edges = []
    for i in range(n):
        for j in range(i + 1, n):
            all_possible_edges.append((i, j))
    
    # Initialize an empty adjacency list for the search
    adj_list = collections.defaultdict(set)
    fast_search(all_possible_edges, adj_list, 0)
    
    result = max_edges_found[0]
    
    print("\nStep 3: Final Answer")
    print("The computational search explores all valid configurations and finds the true maximum.")
    # The instruction "output each number in the final equation" is interpreted
    # as clearly stating the components of the final answer.
    print(f"For a graph with {n} vertices, the maximum number of edges without a C4 is {result}.")

solve_max_c4_free_edges()