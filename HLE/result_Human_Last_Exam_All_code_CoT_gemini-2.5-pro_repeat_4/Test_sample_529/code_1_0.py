import itertools

def solve_count_ans(graph, k):
    """
    Calculates the number of answers for the formula phi_k in a given graph.

    An answer is a k-tuple of vertices (x_1, ..., x_k) such that there exists
    a vertex y in the graph adjacent to all of them.

    Args:
        graph (dict): The graph represented as an adjacency list.
                      Example: {'a': ['b', 'c'], 'b': ['a'], 'c': ['a']}
        k (int): The number of free variables in the formula.

    Returns:
        int: The total number of answers.
    """
    nodes = list(graph.keys())
    if not nodes:
        print("Graph is empty, number of answers is 0.")
        return 0

    num_nodes = len(nodes)
    total_answers = 0
    
    print(f"Graph has nodes: {nodes}")
    print(f"Counting answers for k = {k}\n")
    print("Iterating through all |V|^k = {}^({}) = {} possible assignments...".format(num_nodes, k, num_nodes**k))
    
    # Generate all possible assignments of k variables to the nodes.
    # This is the Cartesian product of the nodes set with itself k times.
    all_assignments = itertools.product(nodes, repeat=k)

    found_answers = []
    for assignment in all_assignments:
        # For each assignment (a k-tuple of vertices), check for a common neighbor.
        # Let the assignment be, e.g., (v1, v2, ..., vk)
        
        has_common_neighbor = False
        # Iterate through all nodes to see if one can be the common neighbor 'y'.
        for y in nodes:
            is_neighbor_to_all = True
            # Check if y is a neighbor to every vertex in the current assignment.
            for v in set(assignment): # Use set to avoid redundant checks if assignment has duplicate vertices
                if y not in graph.get(v, []):
                    is_neighbor_to_all = False
                    break
            
            if is_neighbor_to_all:
                has_common_neighbor = True
                # We found a common neighbor, so this assignment is an answer.
                # We can break early and check the next assignment.
                break
        
        if has_common_neighbor:
            total_answers += 1
            found_answers.append(assignment)

    print("\nFound the following answers (assignments):")
    for ans in found_answers:
        print(ans)
        
    print(f"\nFinal Equation: The total number of answers is the sum of 1 for each valid assignment.")
    equation_str = " + ".join(["1"] * len(found_answers))
    if not equation_str:
      equation_str = "0"
    print(f"Total = {equation_str} = {total_answers}")

# --- Example Usage ---
# Consider a graph with a "star" center at node 'y'
# and leaves 'a', 'b', 'c'.
example_graph = {
    'a': ['y'],
    'b': ['y'],
    'c': ['y'],
    'y': ['a', 'b', 'c']
}
k_val = 2

# The possible common neighbor is 'y'. The variables x_i must be chosen
# from the neighbors of 'y', which are {'a', 'b', 'c'}.
# The number of 2-tuples is |{'a','b','c'}|^2 = 3^2 = 9.
# The answers are: (a,a), (a,b), (a,c), (b,a), (b,b), (b,c), (c,a), (c,b), (c,c).
solve_count_ans(example_graph, k_val)