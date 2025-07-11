import itertools

def get_common_neighbors(graph, vertices):
    """Finds the set of common neighbors for a given set of vertices."""
    if not vertices:
        return set(graph.keys())
    
    # Start with neighbors of the first vertex
    common = set(graph[vertices[0]])
    
    # Intersect with neighbors of the rest of the vertices
    for i in range(1, len(vertices)):
        common.intersection_update(graph[vertices[i]])
        
    return common

def count_answers(graph, k):
    """
    Counts the number of answers for the formula phi_k using the Principle of Inclusion-Exclusion.
    An answer is a k-tuple of vertices (v1, ..., vk) with a common neighbor.
    """
    print(f"Calculating the number of answers for k = {k}")
    print("Formula: Count = |∪ (y in V) N(y)^k|")
    print("Using Principle of Inclusion-Exclusion: Count = Σ (-1)^(|Y|-1) * |N(Y)|^k for non-empty Y ⊆ V\n")

    nodes = sorted(list(graph.keys()))
    n = len(nodes)
    total_count = 0
    
    equation_parts = []

    # Iterate through all non-empty subsets Y of the graph's vertices
    for i in range(1, n + 1):
        term_sum = 0
        for y_subset in itertools.combinations(nodes, i):
            y_subset = list(y_subset)
            
            # Find common neighbors for the current subset Y
            common_neighbors = get_common_neighbors(graph, y_subset)
            num_common = len(common_neighbors)
            
            # Calculate |N(Y)|^k
            term = num_common ** k
            
            # Add or subtract based on the size of Y
            if (i - 1) % 2 == 0: # (-1)^(|Y|-1) is +1
                total_count += term
                term_sum += term
            else: # (-1)^(|Y|-1) is -1
                total_count -= term
                term_sum -= term

        sign = "+" if (i - 1) % 2 == 0 else "-"
        print(f"Term for |Y|={i}: {sign}{abs(term_sum)}")
        if term_sum != 0:
            # Show details for non-zero terms
            for y_subset in itertools.combinations(nodes, i):
                y_subset_str = "{" + ",".join(map(str, y_subset)) + "}"
                common_neighbors = get_common_neighbors(graph, list(y_subset))
                num_common = len(common_neighbors)
                term_val = num_common ** k
                if term_val > 0:
                    print(f"  {sign} |N({y_subset_str})|^k = {sign} {num_common}^{k} = {sign} {term_val}")

    print("\nFinal Equation:")
    print(f"Count = {total_count}")


# Example Graph: A path graph 1-2-3-4
# V = {1, 2, 3, 4}
# E = {(1,2), (2,3), (3,4)}
graph_adj = {
    1: [2],
    2: [1, 3],
    3: [2, 4],
    4: [3]
}

# Parameter k
k_param = 2

count_answers(graph_adj, k_param)