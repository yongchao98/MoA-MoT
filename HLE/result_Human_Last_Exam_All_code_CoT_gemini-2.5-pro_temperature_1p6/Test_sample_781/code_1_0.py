def solve():
    """
    Solves the problem by translating it into a graph theory question.
    The problem is to find the maximum number of edges 'n' in a forest on 'V' vertices
    with no isolated vertices, which is equivalent to finding the number of edges in a tree.
    """
    # The set P contains 5 points, which become the vertices of our graph.
    num_vertices = 5

    # The number of edges in a tree is the number of vertices minus 1.
    # This represents the maximum number of subcontinua in the decomposition.
    n = num_vertices - 1

    # We print the final equation showing the calculation.
    print(f"Let V be the number of special points in the continuum X.")
    print(f"V = {num_vertices}")
    print(f"The problem reduces to finding the maximum number of edges in a tree on V vertices.")
    print(f"The number of edges in a tree is given by the formula n = V - 1.")
    print(f"So, n = {num_vertices} - 1")
    print(f"n = {n}")

solve()