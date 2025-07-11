def solve():
    """
    Calculates the maximum number of colors for a non-complete graph with n vertices.
    """
    # Number of vertices in the graph G
    n = 12345

    # A graph G with n vertices has a chromatic number of n if and only if it is a complete graph (K_n).
    # Since G is not a complete graph, its chromatic number must be strictly less than n.
    # The maximum possible chromatic number is therefore n - 1.
    
    # This maximum is achievable. For instance, the graph K_n with one edge removed
    # is not complete and contains a clique of size n-1, thus its chromatic number is at least n-1.
    # Therefore, the maximum number of colors required is exactly n-1.
    max_colors = n - 1

    print("The number of vertices in the graph is n = 12345.")
    print("For any graph G with n vertices that is not a complete graph, the maximum number of colors needed is n - 1.")
    print("The final calculation is:")
    print(f"{n} - 1 = {max_colors}")

solve()
