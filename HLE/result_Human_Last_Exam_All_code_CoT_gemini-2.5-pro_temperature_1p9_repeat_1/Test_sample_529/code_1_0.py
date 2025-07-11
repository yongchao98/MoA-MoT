import itertools

def count_answers(graph, k):
    """
    Counts the number of answers to the formula phi_k in the given graph using the Principle of Inclusion-Exclusion.

    An answer is a k-tuple of vertices (x_1, ..., x_k) that have a common neighbor.
    Formula: Count = Sum_{W subset V, W!=emptyset} (-1)^(|W|-1) * |N(W)|^k
    where N(W) is the set of common neighbors of all vertices in W.

    Args:
        graph (dict): An adjacency list representation of the graph.
                      Keys are vertices, values are sets of their neighbors.
        k (int): The integer parameter from the problem description.

    Returns:
        int: The number of answers.
    """
    vertices = list(graph.keys())
    n = len(vertices)
    total_answers = 0

    print(f"Graph has vertices: {vertices} and k = {k}")
    print("Calculating using Inclusion-Exclusion: Sum((-1)^{|W|-1} * |N(W)|^{k})")
    print("Non-zero terms of the sum:")
    print("="*30)

    # Iterate over all possible non-empty subset sizes for W
    for i in range(1, n + 1):
        # Iterate over all subsets W of size i
        for W_tuple in itertools.combinations(vertices, i):
            W = set(W_tuple)
            
            # Compute N(W), the set of common neighbors for vertices in W
            if not W:
                continue

            # Start with the neighbors of the first vertex in W
            # and intersect with the neighbors of the rest
            first_vertex = next(iter(W))
            common_neighbors = set(graph.get(first_vertex, set()))
            
            for vertex in W:
                common_neighbors.intersection_update(graph.get(vertex, set()))

            # Calculate the term for the PIE formula
            term = len(common_neighbors) ** k

            # Apply the sign (-1)^(|W|-1)
            if (i - 1) % 2 == 1:
                signed_term = -term
            else:
                signed_term = term

            if signed_term != 0:
                # Per the instructions, outputting each number in the equation.
                # Here we output each non-zero term of the PIE sum.
                sign_char = '+' if signed_term > 0 else '-'
                print(f"W = {W}: {sign_char} |N(W)|^{k} = {sign_char} {len(common_neighbors)}^{k} = {sign_char} {term}")

            total_answers += signed_term
    
    print("="*30)
    print(f"Final Count of Answers: {total_answers}")
    return total_answers

if __name__ == '__main__':
    # Example Graph: A small graph to demonstrate the calculation.
    # 1--2, 1--3, 1--4, 2--3
    # N(1) = {2,3,4}
    # N(2) = {1,3}
    # N(3) = {1,2}
    # N(4) = {1}
    G = {
        1: {2, 3, 4},
        2: {1, 3},
        3: {1, 2},
        4: {1}
    }
    k = 2 # Let's count 2-tuples with a common neighbor
    count_answers(G, k)
