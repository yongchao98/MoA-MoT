import sys

def solve():
    """
    Determines the maximum possible number of different clique sizes
    that can simultaneously appear as induced subgraphs of a single
    graph with n = 128 vertices.
    """
    
    # Number of vertices in the graph
    n = 128

    # The problem asks for the maximum number of different clique sizes that can
    # exist as induced subgraphs in a single graph G with n vertices.
    
    # Step 1: Analyze the nature of induced cliques.
    # If a graph has an induced clique of size k (a K_k), then any subset of
    # j < k vertices from that clique will form an induced clique of size j.
    # This means the set of achievable clique sizes is always {1, 2, ..., k_max},
    # where k_max is the size of the largest induced clique in the graph.
    
    # Step 2: Maximize the number of sizes.
    # To maximize the number of different clique sizes, we must maximize k_max.
    
    # Step 3: Find the maximum possible k_max.
    # For a graph on n vertices, the maximum size of any induced clique is n.
    # This maximum is achieved by the complete graph, K_n.
    
    # Step 4: Construct the solution for n = 128.
    # We choose the graph to be the complete graph on 128 vertices, K_128.
    # In K_128, for any integer k from 1 to 128, any set of k vertices
    # forms an induced clique of size k.
    
    # The maximum number of different clique sizes is therefore n.
    max_sizes = n
    
    # Output the explanation and the final equation.
    print(f"Let the number of vertices be n = {n}.")
    print("To maximize the number of different induced clique sizes, we must maximize the size of the largest induced clique.")
    print(f"The largest possible induced clique in a graph with {n} vertices is of size {n}.")
    print(f"This is achieved by the complete graph K_{n}, which has induced cliques of sizes 1, 2, ..., {n}.")
    print(f"Therefore, the final equation for the maximum number of different clique sizes is: maximum_sizes = {max_sizes}")
    
    # For the final answer extraction format.
    # Suppress all other prints when run in a specific environment for answer extraction.
    if 'ANSWER_ONLY' in sys.argv:
        # Clear previous prints if running in a real environment is complex.
        # For this script, just printing the answer is fine.
        print(f"<<<{max_sizes}>>>")

solve()
