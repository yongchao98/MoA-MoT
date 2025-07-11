import itertools

def find_min_max_models():
    """
    Calculates the number of models for any 3-CNF formula that could have
    generated the given graph via the 3-SAT to Independent Set reduction.
    """
    # Step 1: Define the graph structure from the image.
    # The graph has 12 vertices, which form 4 disjoint triangles (clauses).
    # Let's label the vertices 0 through 11.
    
    # Triangles (Clauses)
    # We can visually partition the vertices into 4 groups of 3.
    # T1 (left): {0, 1, 2}
    # T2 (top-middle): {3, 4, 5}
    # T3 (top-right): {6, 7, 8}
    # T4 (bottom-right): {9, 10, 11}
    clauses = [
        [0, 1, 2],
        [3, 4, 5],
        [6, 7, 8],
        [9, 10, 11]
    ]
    num_clauses = len(clauses)

    # Negation Edges (connections between vertices of different triangles)
    # These are identified by careful visual inspection of the graph image.
    # We use frozenset to represent undirected edges.
    negation_edges = {
        frozenset([2, 3]),
        frozenset([0, 11]),
        frozenset([4, 7]),
        frozenset([5, 6]),
        frozenset([3, 10]),
        frozenset([8, 10]),
        frozenset([7, 9])
    }

    # Step 2: Count the number of independent sets of size k.
    # This corresponds to the number of models for the 3-SAT formula.
    # An independent set of size k=4 must contain exactly one vertex from each clause's triangle.
    
    # Generate all possible selections of one vertex from each triangle.
    all_selections = itertools.product(*clauses)
    
    total_possible_selections = 3**num_clauses
    independent_set_count = 0

    # Iterate through each selection and check if it's an independent set.
    for selection in all_selections:
        is_independent = True
        # Check all pairs of vertices in the current selection for a connecting edge.
        for i in range(num_clauses):
            for j in range(i + 1, num_clauses):
                v1 = selection[i]
                v2 = selection[j]
                if frozenset([v1, v2]) in negation_edges:
                    is_independent = False
                    break
            if not is_independent:
                break
        
        if is_independent:
            independent_set_count += 1

    # Step 3: Determine min and max models.
    # Since the number of models is determined by the fixed graph structure,
    # the minimum and maximum number of models are the same.
    min_models = independent_set_count
    max_models = independent_set_count

    print(f"The graph corresponds to a 3-CNF formula with {num_clauses} clauses.")
    print(f"The total number of possible vertex selections (one from each clause-triangle) is {total_possible_selections}.")
    print(f"After checking for conflicts (negation edges), the number of valid independent sets of size {num_clauses} is {independent_set_count}.")
    print(f"This number is equal to the number of models for any formula that could generate this graph.")
    print(f"Minimum number of models: {min_models}")
    print(f"Maximum number of models: {max_models}")
    print(f"The final answer is the pair (min, max): ({min_models}, {max_models})")

if __name__ == '__main__':
    find_min_max_models()