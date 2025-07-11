import itertools

def count_slices_for_graph(num_vertices, edges, graph_name):
    """
    Counts the number of slices for a given graph.
    This is done by finding the number of valid edge subsets E_A (called colorings)
    and dividing by 2.
    """
    num_edges = len(edges)
    
    valid_colorings = 0
    
    # Iterate through all 2^|E| subsets of edges
    for i in range(1 << num_edges):
        subset_A_edges = []
        for j in range(num_edges):
            if (i >> j) & 1:
                subset_A_edges.append(edges[j])
        
        # An empty edge set is not part of a valid slice partition
        if not subset_A_edges:
            continue

        degrees_A = [0] * num_vertices
        for u, v in subset_A_edges:
            degrees_A[u] += 1
            degrees_A[v] += 1
            
        is_valid = True
        for v_deg in degrees_A:
            # For a cubic graph, the slice condition means the degree in the
            # subgraph must be 1 or 2 (not 0 or 3).
            if not (1 <= v_deg <= 2):
                is_valid = False
                break
        
        if is_valid:
            valid_colorings += 1
            
    num_slices = valid_colorings // 2
    print(f"Graph: {graph_name}, Vertices: {num_vertices}, Edges: {num_edges}")
    print(f"Number of valid edge colorings (C): {valid_colorings}")
    print(f"Number of slices (N): {num_slices}")
    return num_slices

def solve():
    """
    Solves for M(0), M(3), and M(5) based on the plan.
    """
    print("Step 1: Determine M(0)")
    # For N(G) to be a multiple of 0, N(G) must be 0.
    # It is a known theorem that any cubic graph G has N(G) >= 1.
    # Therefore, no such graph exists.
    m0 = "none"
    print(f"Result for M(0): {m0}\n")
    
    print("Step 2: Determine M(3)")
    # Smallest cubic graph is K_4 on 4 vertices.
    # We calculate N(K_4).
    k4_vertices = 4
    k4_edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    n_k4 = count_slices_for_graph(k4_vertices, k4_edges, "K_4")
    
    if n_k4 % 3 == 0:
        m3 = 4
        print(f"N(K_4) = {n_k4} is a multiple of 3. Smallest m is 4.")
    else:
        # This case won't be reached, but is included for completeness.
        m3 = "Could not determine with graphs tested."
    print(f"Result for M(3): {m3}\n")

    print("Step 3: Determine M(5)")
    # Check m=4. N(K_4) = 9 is not a multiple of 5.
    print(f"N(K_4) = {n_k4} is not a multiple of 5. Checking for m=6.")
    # Check m=6. The utility graph K_3,3 is a cubic graph on 6 vertices.
    k33_vertices = 6
    k33_edges = [(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)]
    n_k33 = count_slices_for_graph(k33_vertices, k33_edges, "K_3,3")
    
    if n_k33 % 5 == 0:
        m5 = 6
        print(f"N(K_3,3) = {n_k33} is a multiple of 5. Smallest m is 6.")
    else:
        m5 = "Could not determine with graphs tested."
    print(f"Result for M(5): {m5}\n")
    
    final_answer = f"{m0},{m3},{m5}"
    print("Final Answer format: M(0),M(3),M(5)")
    print(f"Final Answer: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve()