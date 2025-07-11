import itertools

def solve_count_ans(graph, k):
    """
    Calculates the number of answers for the formula phi_k using the
    Principle of Inclusion-Exclusion.

    An answer is a k-tuple of vertices (v1, ..., vk) such that there exists
    a common neighbor y for all vi.

    Args:
        graph (dict): The graph represented as an adjacency list.
                      Example: {'A': {'B'}, 'B': {'A', 'C'}}
        k (int): The number of free variables in the formula.
    """
    vertices = list(graph.keys())
    n = len(vertices)
    total_count = 0
    
    print(f"Calculating the number of answers for k = {k}")
    print("Using the Principle of Inclusion-Exclusion:")
    
    full_equation_terms = []

    for i in range(1, n + 1):
        term_sum = 0
        
        # Get all subsets of vertices of size i
        subsets_Y = itertools.combinations(vertices, i)
        
        subset_terms = []
        for Y in subsets_Y:
            # Find the common neighborhood for the subset Y
            if not Y:
                continue
            
            # Start with the neighbors of the first vertex in Y
            common_neighbors = set(graph[Y[0]])
            
            # Intersect with the neighbors of other vertices in Y
            for j in range(1, len(Y)):
                common_neighbors.intersection_update(graph[Y[j]])
            
            num_common = len(common_neighbors)
            term = num_common ** k
            term_sum += term
            
            # For detailed printing
            subset_terms.append(f"|bigcap(N({','.join(y for y in Y)}))|^{k} = {num_common}^{k} = {term}")


        if i % 2 == 1:
            total_count += term_sum
            sign = "+"
        else:
            total_count -= term_sum
            sign = "-"
        
        if term_sum > 0:
            print(f"Term for |Y|={i}: {sign} ({' + '.join(s.split('= ')[-1] for s in subset_terms if s.endswith(f' = {0}') is False)}) = {sign} {term_sum}")
            full_equation_terms.append(f"{sign} {term_sum}")

    print("\nFinal Equation:")
    print(' '.join(full_equation_terms) + f" = {total_count}")
    print("\nTotal number of answers:")
    print(total_count)


# Example usage:
# A square graph 1-2-3-4-1
example_graph = {
    '1': {'2', '4'},
    '2': {'1', '3'},
    '3': {'2', '4'},
    '4': {'1', '3'}
}
k_val = 2

solve_count_ans(example_graph, k_val)