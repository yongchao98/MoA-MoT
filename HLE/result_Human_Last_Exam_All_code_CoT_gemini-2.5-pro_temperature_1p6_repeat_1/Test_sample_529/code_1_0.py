import itertools

def get_common_neighbors(graph, vertex_subset):
    """
    Finds the set of common neighbors for a subset of vertices.
    """
    if not vertex_subset:
        return set(graph.keys())

    # Get the neighborhood of the first vertex in the subset
    # An explicit copy is needed before starting intersections
    try:
        common_neighbors = set(graph[next(iter(vertex_subset))])
    except StopIteration:
        # This case is handled by the check above, but for safety
        return set(graph.keys())


    # Intersect with the neighborhoods of the other vertices in the subset
    for vertex in vertex_subset:
        if vertex not in graph:
            # If a vertex in the subset does not exist in the graph, it has no neighbors.
            return set()
        common_neighbors.intersection_update(graph[vertex])
    
    return common_neighbors

def powerset(iterable):
    """
    Generates the powerset of an iterable. e.g., powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    """
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s) + 1))

def solve_countans(graph, k):
    """
    Calculates the number of answers for the formula phi_k in graph G.
    It prints the calculation based on the Principle of Inclusion-Exclusion.
    """
    vertices = list(graph.keys())
    total_answers = 0
    
    print(f"Calculating the number of answers for k = {k}")
    print("Using the formula: sum_{Y in P(V), Y non-empty} (-1)^{|Y|-1} * |N(Y)|^k")
    print("-" * 20)

    # Iterate through all non-empty subsets of V
    # Start from subsets of size 1
    for i in range(1, len(vertices) + 1):
        for vertex_subset in itertools.combinations(vertices, i):
            sign = (-1)**(len(vertex_subset) - 1)
            
            common_neighbors = get_common_neighbors(graph, vertex_subset)
            num_common_neighbors = len(common_neighbors)
            
            term_value = sign * (num_common_neighbors ** k)
            total_answers += term_value
            
            sign_char = '+' if sign == 1 else '-'
            print(f"Y = {set(vertex_subset)}: {sign_char} |N(Y)|^k = {sign_char} ({num_common_neighbors})^{k} = {term_value}")

    print("-" * 20)
    print(f"The total number of answers is: {total_answers}")

if __name__ == '__main__':
    # Example Graph G: A path of 4 vertices: 1 -- 2 -- 3 -- 4
    # The graph is represented as an adjacency list (dictionary in Python)
    G = {
        1: [2],
        2: [1, 3],
        3: [2, 4],
        4: [3]
    }
    k_param = 2
    
    solve_countans(G, k_param)

<<<B>>>