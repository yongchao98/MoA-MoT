def solve_max_clique_sizes():
    """
    This function determines the maximum possible number of different clique sizes
    that can simultaneously appear as induced subgraphs of a single graph with n vertices.

    The logic is as follows:
    1.  The set of induced clique sizes in any graph G is {1, 2, ..., omega(G)},
        where omega(G) is the size of the largest clique in G.
    2.  Thus, the number of different clique sizes is omega(G).
    3.  The problem is to find the maximum possible value of omega(G) for a graph
        on n vertices.
    4.  The clique size is limited by the number of vertices, so omega(G) <= n.
    5.  The complete graph K_n has a clique of size n, so the maximum is n.
    """

    # The number of vertices in the graph is given.
    n = 128

    # The maximum possible number of different clique sizes is n.
    max_sizes = n

    # The final equation is simply max_sizes = n.
    # The number we are interested in is the result.
    print(max_sizes)

solve_max_clique_sizes()