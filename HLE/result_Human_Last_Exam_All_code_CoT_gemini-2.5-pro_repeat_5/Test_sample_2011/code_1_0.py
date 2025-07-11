import sys

def solve_clique_sizes():
    """
    Determines the maximum possible number of different clique sizes
    that can simultaneously appear as induced subgraphs of a single graph
    on n vertices.
    """
    # The number of vertices in the graph, as specified by the problem.
    n = 128

    print("Step 1: Understand the problem.")
    print("We want to find the maximum number of different values 'k' such that a clique of size k (K_k)")
    print(f"can be found as an induced subgraph in a single graph G with n = {n} vertices.\n")

    print("Step 2: Analyze the property of induced cliques.")
    print("Let's assume a graph G has an induced clique of size k_max. This means there is a set of k_max vertices where every vertex is connected to every other vertex in the set.")
    print("If we take any subset of k vertices from this set (where 1 <= k <= k_max), this subset also forms a complete graph.")
    print("This means that if K_{k_max} is an induced subgraph, then K_k is also an induced subgraph for all k in {1, 2, ..., k_max}.")
    print("Therefore, the set of all possible induced clique sizes in any graph is always of the form {1, 2, ..., k_max}, where k_max is the size of the largest induced clique in the graph.\n")

    print("Step 3: Simplify the problem.")
    print("The number of different clique sizes is simply k_max.")
    print(f"The problem is now to find the maximum possible value of k_max for a graph on n = {n} vertices.\n")

    print("Step 4: Determine the theoretical maximum for k_max.")
    print("A clique is a set of vertices from the graph. The size of a clique cannot be larger than the total number of vertices.")
    print(f"Therefore, k_max must be less than or equal to n. So, k_max <= {n}.\n")

    print("Step 5: Show that this maximum is achievable.")
    print(f"Consider the complete graph K_n (in our case, K_{n}). This is a graph with n vertices where every vertex is connected to every other vertex.")
    print(f"In K_{n}, any subset of k vertices forms an induced clique of size k.")
    print(f"This is true for all k from 1 to n. The largest clique is the graph itself, with size n.")
    print(f"So, for the graph K_{n}, k_max = n is achieved.\n")

    print("Step 6: Conclude and state the final answer.")
    print("The maximum possible number of different clique sizes is the maximum possible value of k_max, which is n.")
    
    # Final calculation
    max_sizes = n
    
    # Output the final "equation" as requested
    print("\n--- Final Answer ---")
    print(f"Let m be the maximum number of different clique sizes.")
    print(f"Let n be the number of vertices.")
    print(f"The final equation is: m = n")
    print(f"For n = {n}, the result is:")
    print(f"m = {max_sizes}")


solve_clique_sizes()

# Return the final numerical answer in the specified format
sys.stdout.write("\n<<<128>>>\n")