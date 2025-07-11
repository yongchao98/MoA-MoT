def solve_max_clique_sizes():
    """
    Determines and explains the maximum possible number of different clique sizes
    in a graph on n vertices.
    """
    # The number of vertices in the graph.
    n = 128

    # The maximum possible number of different clique sizes is determined by logic.
    # 1. The size 'k' of any clique in a graph on 'n' vertices cannot exceed 'n'.
    #    Therefore, the number of different clique sizes is at most 'n'.
    # 2. The complete graph K_n is a graph on 'n' vertices.
    # 3. In K_n, any subset of 'k' vertices (for 1 <= k <= n) induces a clique of size 'k'.
    # 4. This means K_n contains induced cliques for all sizes from 1 to n.
    # 5. Since the maximum is at most 'n' and we found a graph that achieves 'n', the maximum is 'n'.

    max_sizes = n

    # The final equation is simply that the maximum number equals n.
    print(f"The number of vertices is n = {n}.")
    print(f"The maximum possible number of different clique sizes in a graph on {n} vertices is {max_sizes}.")
    print(f"This is achieved using a complete graph K_{n}, which contains induced cliques of sizes 1, 2, ..., {n}.")

solve_max_clique_sizes()