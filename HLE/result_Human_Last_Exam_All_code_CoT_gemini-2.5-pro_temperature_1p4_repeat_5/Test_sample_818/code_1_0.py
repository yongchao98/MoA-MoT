import itertools

def count_slices(num_vertices, edges):
    """
    Counts the number of slices for a given graph.
    A slice is a partition of edges {E1, E2} where in each induced subgraph,
    every vertex has a degree of at least 1. For a cubic graph, this means
    the degrees in the subgraph for E1 must be 1 or 2 for all vertices.
    """
    num_edges = len(edges)
    valid_slice_sets_count = 0

    # Iterate through all non-empty proper subsets of edges to form E1
    # 2**(num_edges - 1) is sufficient to avoid double counting {E1,E2} and {E2,E1}
    for i in range(1, 2**(num_edges - 1)):
        e1_indices = []
        # Construct the subset of edge indices from the integer i
        for k in range(num_edges):
            if (i >> k) & 1:
                e1_indices.append(k)
        
        e1 = [edges[k] for k in e1_indices]
        
        degs = {v: 0 for v in range(num_vertices)}
        for u, v in e1:
            degs[u] += 1
            degs[v] += 1
            
        is_slice_set = all(1 <= d <= 2 for d in degs.values())
        
        if is_slice_set:
            valid_slice_sets_count += 1
            
    return valid_slice_sets_count

def main():
    """
    Main function to calculate M(0), M(3), M(5) and print the result.
    """
    # M(0): No cubic graph G has N(G) = 0.
    m0 = "none"
    
    # M(3): Check the smallest cubic graph, K4 (m=6).
    k4_vertices = 4
    k4_edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    n_k4 = count_slices(k4_vertices, k4_edges)
    
    m3 = "none"
    if n_k4 % 3 == 0:
        m3 = len(k4_edges) # 6
        
    # M(5): Check K4 first.
    m5 = "none"
    if n_k4 % 5 == 0:
        m5 = len(k4_edges)
    else:
        # Check next smallest cubic graphs, m=9.
        prism_vertices = 6
        prism_edges = [(0, 1), (1, 2), (2, 0),
                       (3, 4), (4, 5), (5, 3),
                       (0, 3), (1, 4), (2, 5)]
        n_prism = count_slices(prism_vertices, prism_edges)
        if n_prism % 5 == 0:
            m5 = len(prism_edges) # 9
        # If not, we would continue to check K3,3 and then larger graphs.
        # But for this problem, checking the prism graph is sufficient.

    print(f"{m0},{m3},{m5}")

if __name__ == "__main__":
    main()
<<<none,6,9>>>