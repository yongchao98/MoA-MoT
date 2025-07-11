def solve():
    """
    Calculates the numbers of isomorphism classes of vertex-transitive graphs
    with 8 vertices and vertex degree j for j=0, ..., 7.
    """
    # These values are based on known classifications from graph theory.
    # n_j is the number of non-isomorphic vertex-transitive graphs of degree j on 8 vertices.
    
    # n_0: The empty graph (8K_1)
    n_0 = 1
    
    # n_1: A perfect matching (4K_2)
    n_1 = 1
    
    # n_2: The cycle C_8 and the disjoint union 2C_4
    n_2 = 2
    
    # n_3: The cube graph Q_3, the circulant graph C_8(1,4), and 2K_4
    n_3 = 3
    
    # For a graph G of degree j on n vertices, its complement has degree n-1-j.
    # For n=8, this means n_j = n_{7-j}.
    n_4 = n_3
    n_5 = n_2
    n_6 = n_1
    n_7 = n_0
    
    result = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]
    print(f"{result}")

solve()