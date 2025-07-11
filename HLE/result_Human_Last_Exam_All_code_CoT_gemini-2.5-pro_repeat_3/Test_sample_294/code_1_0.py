import itertools

def count_k_matchings(edges, k):
    """
    Counts the number of k-matchings in a graph using a simple brute-force approach.

    A k-matching is a set of k edges with no shared vertices. This function
    checks every combination of k edges to see if it forms a valid matching.

    Note: This implementation is for demonstration and is not time-efficient.
    Its complexity is O(|E|^k), which is not subcubic for k>=2.
    The theoretical possibility of subcubic algorithms for k=3 relies on
    much more advanced techniques like fast matrix multiplication.

    Args:
        edges (list of tuples): The list of edges in the graph.
        k (int): The size of the matching to count.

    Returns:
        int: The total number of k-matchings.
    """
    if k == 0:
        return 1
    
    # Get all combinations of k edges from the edge set
    edge_combinations = itertools.combinations(edges, k)
    
    matching_count = 0
    for combo in edge_combinations:
        # A combo is a potential k-matching, e.g., ((1,2), (3,4), (5,6)) for k=3
        
        # To check if it's a valid matching, we verify that no vertices are shared.
        # We can do this by flattening the list of vertices and checking for uniqueness.
        vertices_in_combo = [vertex for edge in combo for vertex in edge]
        
        # If the number of vertices is equal to the number of unique vertices,
        # it means there are no duplicates, and it's a valid matching.
        # A k-matching must involve 2*k distinct vertices.
        if len(set(vertices_in_combo)) == 2 * k:
            matching_count += 1
            
    return matching_count

if __name__ == '__main__':
    # The question is: what is the maximum k for which k-matchings can be
    # counted in subcubic time (O(|V|^(3-eps)))?
    
    # Based on fine-grained complexity theory:
    # - For k=1, 2, 3: YES, subcubic algorithms exist.
    #   The algorithm for k=3 is complex and related to triangle counting (O(|V|^omega)).
    # - For k=4: NO, it is assumed not to be possible in subcubic time, as it's at
    #   least as hard as counting 4-cliques.
    
    # Therefore, the maximum value is 3.
    max_k = 3

    print(f"The final answer is that the maximum value of k is {max_k}.")
    print("---")

    # Let's use the code to count the number of 3-matchings in an example graph.
    # Consider a complete graph on 6 vertices (K6). A 3-matching would partition the vertices.
    num_vertices_example = 6
    vertices = range(1, num_vertices_example + 1)
    
    # Generate all edges for a complete graph K6
    example_edges = list(itertools.combinations(vertices, 2))
    
    # Count the number of 3-matchings in K6
    num_matchings = count_k_matchings(example_edges, max_k)
    
    print(f"As an example, we count the number of {max_k}-matchings in a complete graph with {num_vertices_example} vertices.")
    print(f"The number of {max_k}-matchings is: {num_matchings}")
    
    # The equation for the number of perfect matchings in K_{2n} is (2n-1)!!
    # For n=3 (K6), this is 5!! = 5 * 3 * 1 = 15. Our code should find this.
    print(f"The final equation for a perfect matching in K_6 is (2*{max_k}-1)!! = 5*3*1 = 15.")
