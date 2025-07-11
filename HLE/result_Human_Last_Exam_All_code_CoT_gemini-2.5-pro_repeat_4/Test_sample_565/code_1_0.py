def solve():
    """
    Calculates the numbers n_j of isomorphism classes of vertex-transitive graphs
    with 8 vertices and vertex degree precisely j for j=0, ..., 7.
    """
    # n_j is the number of non-isomorphic vertex-transitive graphs of degree j.
    # For n=8 vertices, n_j = n_{7-j}.

    # n_0: The empty graph on 8 vertices. It is unique and vertex-transitive.
    n_0 = 1
    
    # n_1: A 1-regular graph on 8 vertices is a perfect matching (4K_2).
    # It is unique up to isomorphism and vertex-transitive.
    n_1 = 1
    
    # n_2: A 2-regular graph is a union of cycles.
    # On 8 vertices, only C_8 and 2C_4 are vertex-transitive.
    n_2 = 2
    
    # n_3: There are 4 cubic vertex-transitive graphs on 8 vertices.
    # Three are connected (the Cube graph, the Mobius ladder M_8, and one other)
    # and one is disconnected (2K_4).
    n_3 = 4
    
    # Using symmetry n_j = n_{7-j}
    n_4 = n_3
    n_5 = n_2
    n_6 = n_1
    n_7 = n_0
    
    result = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]
    
    print(result)

solve()