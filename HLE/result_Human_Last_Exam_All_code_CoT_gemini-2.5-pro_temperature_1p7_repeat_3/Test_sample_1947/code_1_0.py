def solve_coefficients():
    """
    This function determines the coefficients c_1, c_2, c_3, c_4, c_5 based on combinatorial analysis.

    The total number of closed tree-like walks of length 6 in a graph X is the sum of
    such walks over all of its distinct tree subgraphs. The possible underlying trees
    must have 1, 2, or 3 edges.

    Let N(T) be the number of distinct closed tree-like walks of length 6 whose edge set
    forms a tree isomorphic to T.
    Through combinatorial counting, we find:
    - N(P_2) = 2  (P_2 is a single edge)
    - N(P_3) = 12 (P_3 is a path with 2 edges)
    - N(P_4) = 6  (P_4 is a path with 3 edges)
    - N(K_1,3) = 12 (K_1,3 is a star graph with 3 edges)

    The total count in a graph X is given by the sum over its subgraphs:
    Total_Walks = N(P_2)*n(P_2) + N(P_3)*n(P_3) + N(P_4)*n(P_4) + N(K_1,3)*n(K_1,3)
    where n(T) is the number of subgraphs in X isomorphic to T.

    These subgraph counts are related to the terms in the given expression:
    - n(P_2) = e
    - n(P_3) = sum(deg(v) choose 2) over all vertices v
    - n(P_4) = p
    - n(K_1,3) = sum(deg(v) choose 3) over all vertices v

    Substituting these into the total walks formula gives:
    Total_Walks = 2*e + 12*sum(deg C 2) + 6*p + 12*sum(deg C 3)

    Comparing this with the given expression:
    c_1*e + c_2*k + c_3*p + c_4*sum(deg C 2) + c_5*sum(deg C 3)
    allows us to determine the coefficients by matching terms.
    """
    
    # c_1 is the coefficient for e (number of edges, i.e., P_2 subgraphs)
    c1 = 2
    
    # c_2 is the coefficient for k (number of K_3 subgraphs)
    # Tree-like walks do not depend on non-tree structures like triangles in a way
    # that isn't already captured by their tree subgraphs. So, the coefficient is 0.
    c2 = 0
    
    # c_3 is the coefficient for p (number of P_4 subgraphs)
    c3 = 6
    
    # c_4 is the coefficient for sum(deg(v) choose 2) (number of P_3 subgraphs)
    c4 = 12
    
    # c_5 is the coefficient for sum(deg(v) choose 3) (number of K_1,3 subgraphs)
    c5 = 12
    
    print(c1)
    print(c2)
    print(c3)
    print(c4)
    print(c5)

solve_coefficients()