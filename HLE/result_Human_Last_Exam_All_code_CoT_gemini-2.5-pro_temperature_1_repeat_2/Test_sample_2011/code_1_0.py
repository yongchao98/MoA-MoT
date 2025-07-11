import sys

def solve_clique_sizes():
    """
    Determines the maximum possible number of different clique sizes
    in a graph on n vertices.
    """
    # The number of vertices in the graph.
    n = 128

    # Step 1: Understand the problem and establish an upper bound.
    # A "clique" in a graph is a subset of its vertices where every two distinct
    # vertices are adjacent. A clique of size k is equivalent to an induced
    # subgraph isomorphic to the complete graph K_k.
    # The question asks for the maximum number of different values of k
    # for which a K_k can be found as an induced subgraph in a single graph G
    # with n=128 vertices.

    # The size of any clique, k, must be a positive integer.
    # Also, a clique is a subset of the graph's vertices, so its size cannot
    # be larger than the total number of vertices, n.
    # Therefore, any possible clique size k must satisfy 1 <= k <= n.
    # This means the set of all possible clique sizes is a subset of {1, 2, ..., n}.
    # Consequently, the number of different clique sizes cannot exceed n.
    upper_bound = n

    # Step 2: Show that this upper bound is achievable.
    # To do this, we need to construct a graph on n vertices that contains
    # cliques of all sizes from 1 to n.
    # Consider the complete graph on n vertices, K_n. In this graph, every
    # pair of distinct vertices is connected by an edge.

    # Let's verify the clique sizes in K_n.
    # For any integer k where 1 <= k <= n, we can choose any subset of k vertices
    # from the n vertices of K_n.
    # The subgraph induced by this subset of k vertices will include all edges
    # between them, because all possible edges exist in K_n.
    # Thus, this induced subgraph is a complete graph of size k, or K_k.

    # This proves that the graph K_n contains cliques of all sizes k = 1, 2, ..., n.
    # The number of different clique sizes in K_n is therefore n.

    # Step 3: Conclude the result.
    # Since the number of different clique sizes cannot exceed n, and we have
    # found a graph (K_n) that achieves exactly n different clique sizes, the
    # maximum possible number is n.

    # The final "equation" is simply that the maximum number of sizes is equal
    # to the number of vertices.
    max_sizes = n
    
    print(f"The number of vertices is n = {n}.")
    print(f"The maximum possible number of different clique sizes in a graph on n vertices is n.")
    print("This is because the complete graph K_n contains induced subgraphs of every size from 1 to n.")
    print("\nFinal Equation:")
    print(f"max_possible_sizes = n")
    print(f"max_possible_sizes = {max_sizes}")

solve_clique_sizes()