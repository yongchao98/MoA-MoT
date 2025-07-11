def solve_coloring_problem():
    """
    Calculates the maximum number of colours needed for a graph with n vertices
    that is not a complete graph.
    """
    # Number of vertices in the graph G
    n = 12345

    # According to Brooks's Theorem, for a graph G that is not a complete graph
    # or an odd cycle, the chromatic number chi(G) is at most the maximum degree Delta(G).
    # chi(G) <= Delta(G).

    # For any graph with n vertices, the maximum degree is at most n-1.
    # So, chi(G) <= Delta(G) <= n-1.

    # An odd cycle graph with n vertices would require 3 colours. Since n=12345 is odd,
    # this is a possibility, yielding 3.

    # We need to find the maximum possible number of colors. To do this, we try to
    # maximize chi(G). The upper bound is n-1. We can check if this bound is achievable.

    # Consider a graph G constructed by taking a complete graph K_n and removing
    # one edge. This graph is not complete. It has a clique of size n-1, so it
    # requires at least n-1 colours. Since chi(G) <= n-1, its chromatic number is
    # exactly n-1.

    # The maximum of the possible values (3 for an odd cycle, and n-1 for the constructed graph)
    # gives the answer.
    max_colors = n - 1

    # Print the final equation as requested.
    print(f"The maximum number of colors needed is given by the equation:")
    print(f"{n} - 1 = {max_colors}")

solve_coloring_problem()