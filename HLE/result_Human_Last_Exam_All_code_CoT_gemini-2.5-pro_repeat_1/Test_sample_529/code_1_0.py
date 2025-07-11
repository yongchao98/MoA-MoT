import itertools

def solve_count_ans(graph, k):
    """
    Calculates the number of answers for the formula phi_k using the Principle of Inclusion-Exclusion.

    An answer is a k-tuple of vertices (v_1, ..., v_k) that have a common neighbor.

    Args:
        graph (dict): The graph represented as an adjacency list.
                      Keys are vertices (int), values are sets of neighbors.
        k (int): The parameter k for the formula phi_k.
    """
    vertices = list(graph.keys())
    num_vertices = len(vertices)
    total_count = 0

    # The equation is: Count = sum_{S subset V, S non-empty} (-1)^{|S|-1} * |intersection_{y in S} N(y)|^k
    # We will iterate through all possible sizes of subsets S
    for subset_size in range(1, num_vertices + 1):
        # Contribution for subsets of the current size
        term_sum = 0
        
        # Iterate through all subsets S of the current size
        for s_indices in itertools.combinations(range(num_vertices), subset_size):
            subset_S = [vertices[i] for i in s_indices]

            # Calculate the intersection of neighborhoods for vertices in subset_S
            if not subset_S:
                continue

            # Start with the neighborhood of the first vertex in S
            common_neighbors = set(graph[subset_S[0]])
            
            # Intersect with the neighborhoods of the other vertices in S
            for i in range(1, len(subset_S)):
                vertex = subset_S[i]
                common_neighbors.intersection_update(graph[vertex])
            
            num_common_neighbors = len(common_neighbors)
            term = pow(num_common_neighbors, k)
            term_sum += term
            
            # For demonstration, print the calculation for this term of the equation
            intersection_str = " & ".join([f"N({v})" for v in subset_S])
            print(f"Term for S={set(subset_S)}: |{intersection_str}| = {num_common_neighbors}. Contribution: {num_common_neighbors}^{k} = {term}")


        # Apply the sign from the PIE formula: (-1)^{|S|-1}
        sign = (-1)**(subset_size - 1)
        
        if sign == 1:
            print(f"\nAdding contributions for subsets of size {subset_size}: {term_sum}")
            total_count += term_sum
        else:
            print(f"\nSubtracting contributions for subsets of size {subset_size}: {term_sum}")
            total_count -= term_sum
        
        print(f"Current total count = {total_count}\n")


    print("="*20)
    print(f"Final number of answers for phi_{k} is: {total_count}")

# Example Usage:
if __name__ == '__main__':
    # A simple path graph: 0 -- 1 -- 2
    # V = {0, 1, 2}
    # E = {(0,1), (1,0), (1,2), (2,1)}
    # N(0) = {1}, N(1) = {0, 2}, N(2) = {1}
    example_graph = {
        0: {1},
        1: {0, 2},
        2: {1}
    }
    k_param = 2

    # For G={0-1-2}, k=2:
    # Answers with witness y=0: N(0)={1}. Tuples: (1,1). Set={(1,1)}
    # Answers with witness y=1: N(1)={0,2}. Tuples: (0,0),(0,2),(2,0),(2,2). Set={(0,0),(0,2),(2,0),(2,2)}
    # Answers with witness y=2: N(2)={1}. Tuples: (1,1). Set={(1,1)}
    # Union of these sets: {(1,1), (0,0), (0,2), (2,0), (2,2)}. Total size = 5.
    
    solve_count_ans(example_graph, k_param)