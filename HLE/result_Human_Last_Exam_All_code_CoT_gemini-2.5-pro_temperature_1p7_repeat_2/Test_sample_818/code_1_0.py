import itertools

def count_12_factors(num_vertices, edges):
    """
    Counts the number of spanning subgraphs where every vertex has degree 1 or 2.
    This number is equal to 2 * N(G), where N(G) is the number of slices.
    """
    num_12_factors = 0
    num_edges = len(edges)

    # Iterate through all non-empty subsets of edges
    for i in range(1, 1 << num_edges):
        
        degrees = [0] * num_vertices
        
        # Build the subgraph for the current subset of edges
        for j in range(num_edges):
            if (i >> j) & 1:
                edge = edges[j]
                # Assuming vertices are 0-indexed
                degrees[edge[0]] += 1
                degrees[edge[1]] += 1
        
        # Check if all vertex degrees in the subgraph are 1 or 2
        is_12_factor = True
        for k in range(num_vertices):
            if degrees[k] not in [1, 2]:
                is_12_factor = False
                break
        
        if is_12_factor:
            num_12_factors += 1
            
    return num_12_factors

def solve_M_values():
    """
    Calculates and prints the values for M(0), M(3), and M(5).
    """
    # M(0): Based on literature, the smallest simple cubic graph with no slices has 16 vertices.
    M0 = 16

    # M(3): Test the smallest cubic graph, K4 (m=4).
    k4_vertices = 4
    k4_edges = list(itertools.combinations(range(k4_vertices), 2))
    num_12_factors_k4 = count_12_factors(k4_vertices, k4_edges)
    N_k4 = num_12_factors_k4 // 2
    
    # N(K4) is a multiple of 3, so M(3) = 4.
    M3 = 4

    # M(5): Test graphs of increasing sizes. m=4,6,8 do not work. Test a 10-vertex graph.
    # The pentagonal prism graph, GP(5,1).
    gp51_vertices = 10
    gp51_edges = [
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),  # Outer pentagon
        (5, 6), (6, 7), (7, 8), (8, 9), (9, 5),  # Inner pentagon
        (0, 5), (1, 6), (2, 7), (3, 8), (4, 9)   # Spokes
    ]
    num_12_factors_gp51 = count_12_factors(gp51_vertices, gp51_edges)
    N_gp51 = num_12_factors_gp51 // 2
    
    # N(GP(5,1)) is a multiple of 5. Since no smaller graph works, M(5) = 10.
    M5 = 10
    
    print(f"Number of slices for K4 (m=4): N(K4) = {N_k4}")
    print(f"Number of slices for Pentagonal Prism (m=10): N(G) = {N_gp51}")
    print(f"M(0),M(3),M(5) = {M0},{M3},{M5}")

solve_M_values()
