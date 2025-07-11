def solve_graph_count():
    """
    Calculates the numbers n_j of isomorphism classes of vertex-transitive
    graphs with 8 vertices and vertex degree j for j=0, ..., 7.
    """
    # Based on the known classification of vertex-transitive graphs on 8 vertices.
    
    # n_0: The empty graph (8K1).
    n_0 = 1
    
    # n_1: The perfect matching graph (4K2).
    n_1 = 1
    
    # n_2: The cycle graph C8 and two disjoint C4 graphs (2C4).
    n_2 = 2
    
    # n_3: The Cubical graph and the circulant graph C(8,{1,4}).
    n_3 = 2
    
    # The number of vertex-transitive graphs of degree j on n vertices is the
    # same as the number of those with degree n-1-j.
    # For n=8, n_j = n_{7-j}.
    
    # n_4 is the same as n_3
    n_4 = n_3
    
    # n_5 is the same as n_2
    n_5 = n_2
    
    # n_6 is the same as n_1
    n_6 = n_1
    
    # n_7 is the same as n_0 (the complete graph K8).
    n_7 = n_0
    
    result = [n_0, n_1, n_2, n_3, n_4, n_5, n_6, n_7]
    
    print(f"The number of isomorphism classes for degree j=0 is n_0 = {n_0}")
    print(f"The number of isomorphism classes for degree j=1 is n_1 = {n_1}")
    print(f"The number of isomorphism classes for degree j=2 is n_2 = {n_2}")
    print(f"The number of isomorphism classes for degree j=3 is n_3 = {n_3}")
    print(f"The number of isomorphism classes for degree j=4 is n_4 = {n_4}")
    print(f"The number of isomorphism classes for degree j=5 is n_5 = {n_5}")
    print(f"The number of isomorphism classes for degree j=6 is n_6 = {n_6}")
    print(f"The number of isomorphism classes for degree j=7 is n_7 = {n_7}")
    
    print("\nFinal list in the format [n_0, n_1, ..., n_7]:")
    print(result)

solve_graph_count()