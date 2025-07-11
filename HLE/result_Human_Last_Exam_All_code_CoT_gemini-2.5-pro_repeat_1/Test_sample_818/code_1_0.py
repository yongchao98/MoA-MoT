import itertools

def calculate_n_g(num_vertices, edges):
    """
    Calculates N(G) for a given graph G.
    N(G) is the number of "slices". A slice is a partition of edges into two
    classes where each vertex is incident to at least one edge in each class.
    
    This is done by checking all 2^E edge colorings, counting the valid ones
    (where no vertex is monochromatic), and dividing by 2.
    """
    num_edges = len(edges)
    
    # Create an incidence map for easy lookup: vertex -> list of edge indices
    incidence_map = {v: [] for v in range(num_vertices)}
    for i, edge in enumerate(edges):
        incidence_map[edge[0]].append(i)
        incidence_map[edge[1]].append(i)

    valid_s_colorings = 0
    # Iterate through all 2^num_edges possible 2-colorings of the edges
    # A bit '0' can be class A, a bit '1' can be class B
    for i in range(2**num_edges):
        is_valid_coloring = True
        # Check each vertex for the slice condition
        for v in range(num_vertices):
            edge_indices = incidence_map[v]
            # Skip isolated vertices, although none exist in cubic graphs
            if not edge_indices:
                continue

            # Get the color of the first incident edge
            first_edge_color = (i >> edge_indices[0]) & 1
            
            # Check if all other incident edges have the same color
            is_monochromatic = True
            for edge_idx in edge_indices[1:]:
                if (i >> edge_idx) & 1 != first_edge_color:
                    is_monochromatic = False
                    break
            
            if is_monochromatic:
                is_valid_coloring = False
                break  # This coloring is invalid, move to the next one
        
        if is_valid_coloring:
            valid_s_colorings += 1
            
    # N(G) is the number of partitions. Each valid coloring C and its
    # complement C' form a single partition, so we divide by 2.
    return valid_s_colorings // 2

def solve_and_print():
    """
    Solves for M(0), M(3), and M(5) and prints the results and reasoning.
    """
    # Part 1: Determine M(0)
    # M(0) is the smallest m for which N(G) is a multiple of 0, meaning N(G)=0.
    # The existence of a cubic graph G with N(G)=0 is an open problem.
    # Based on the conjecture that N(G) > 0 for all cubic graphs, M(0) is "none".
    m0_val = "none"
    print("Step 1: Determine M(0)")
    print("M(0) requires a cubic graph G with N(G) = 0.")
    print("The existence of such a graph is an open mathematical problem. It is conjectured that none exists.")
    print(f"Therefore, M(0) is taken to be {m0_val}.\n")

    # Part 2: Determine M(3)
    # Search for the smallest m where N(G) is a multiple of 3.
    # The smallest cubic graph has m=4 (K_4).
    print("Step 2: Determine M(3)")
    print("Searching for the smallest m where N(G) is a multiple of 3.")
    
    v_k4 = 4
    e_k4 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    n_k4 = calculate_n_g(v_k4, e_k4)
    
    print(f"For m = 4, the graph is K_4. The number of slices N(K_4) is {n_k4}.")
    
    m3_val = "none"
    if n_k4 % 3 == 0:
        m3_val = 4
        print(f"The equation is N(K_4) = {n_k4}. Since {n_k4} is a multiple of 3, and m=4 is the smallest possible size for a cubic graph, M(3) = {m3_val}.\n")
    
    # Part 3: Determine M(5)
    # Search for the smallest m where N(G) is a multiple of 5.
    print("Step 3: Determine M(5)")
    print("Searching for the smallest m where N(G) is a multiple of 5.")
    
    m5_val = "none"
    print(f"For m = 4 (K_4), we found N(K_4) = {n_k4}.")
    if n_k4 % 5 == 0:
        m5_val = 4
        print(f"N(K_4) is a multiple of 5, so M(5) = 4.")
    else:
        print(f"N(K_4) = {n_k4} is not a multiple of 5. We check the next possible size, m=6.")
        
        # Test K_3,3 (utility graph), one of the two cubic graphs on 6 vertices.
        v_k33 = 6
        e_k33 = [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)]
        n_k33 = calculate_n_g(v_k33, e_k33)
        print(f"For m = 6, testing the graph K_3,3. The number of slices N(K_3,3) is {n_k33}.")

        if n_k33 % 5 == 0:
            m5_val = 6
            print(f"The equation is N(K_3,3) = {n_k33}. Since {n_k33} is a multiple of 5, and m=4 did not work, M(5) = {m5_val}.\n")
        else:
             print(f"N(K_3,3) = {n_k33} is not a multiple of 5. We would test the other m=6 graph if necessary.")

    # Final Answer
    final_answer = f"{m0_val},{m3_val},{m5_val}"
    print(f"<<<{final_answer}>>>")

solve_and_print()