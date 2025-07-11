import itertools

def solve_count_ans(graph, k):
    """
    Solves the CountAns problem for a given graph and parameter k using a brute-force approach.

    Args:
        graph (dict): An adjacency list representation of the graph.
                      Example: {0: [1, 2], 1: [0], 2: [0]}
        k (int): The parameter k for the formula phi_k.
    """
    vertices = list(graph.keys())
    n = len(vertices)
    if n == 0:
        print("Graph is empty. Number of answers is 0.")
        return

    # Generate all k-tuples of vertices with replacement
    all_k_tuples = itertools.product(vertices, repeat=k)
    
    answer_count = 0
    answer_terms = []

    for tpl in all_k_tuples:
        # For each tuple, check for a common neighbor
        has_common_neighbor = False
        for y in vertices:
            # Check if y is a neighbor to all vertices in the tuple
            is_common_neighbor = True
            # Use a set for faster lookup of vertices in the tuple
            v_set = set(tpl)
            for v in v_set:
                if y not in graph.get(v, []):
                    is_common_neighbor = False
                    break
            
            if is_common_neighbor:
                has_common_neighbor = True
                break # Found a common neighbor, no need to check others for this tuple
        
        if has_common_neighbor:
            # print(f"Tuple {tpl} is an answer.")
            answer_count += 1
            answer_terms.append("1")

    if not answer_terms:
        print("There are no answers for the given graph and k.")
        print("Final equation: 0 = 0")
    else:
        equation_str = " + ".join(answer_terms)
        print(f"Found {answer_count} answers.")
        print(f"The final equation is: {equation_str} = {answer_count}")


# Example Usage:
# Consider a graph where vertices {0, 1} have common neighbors {2, 3}.
# G = (V, E) with V = {0, 1, 2, 3} and E = {(0,2), (2,0), (0,3), (3,0), (1,2), (2,1), (1,3), (3,1)}
example_graph = {
    0: [2, 3],
    1: [2, 3],
    2: [0, 1],
    3: [0, 1]
}
k_param = 2

# For k=2, the answers are tuples (v1, v2) where v1, v2 are from {0, 1}.
# These are (0,0), (0,1), (1,0), (1,1). There should be 4 answers.
solve_count_ans(example_graph, k_param)
