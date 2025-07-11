def solve():
    """
    This function calculates and prints the coefficients for the number of closed tree-like walks of length 6.
    The expression is of the form:
    c_1 * e + c_2 * k + c_3 * p + c_4 * sum(deg(v) choose 2) + c_5 * sum(deg(v) choose 3)
    where e is number of edges, k is number of K_3, p is number of P_4.
    """
    
    # Based on the combinatorial analysis of tree-like walks:
    # c1 corresponds to walks on a P_2 tree structure.
    c1 = 2
    
    # c2 corresponds to walks involving a K_3. By definition, tree-like walks cannot
    # have a cycle in their edge set, so this term's coefficient is 0.
    c2 = 0
    
    # c3 corresponds to walks on a P_4 tree structure.
    c3 = 6
    
    # c4 corresponds to walks on a P_3 tree structure.
    c4 = 12
    
    # c5 corresponds to walks on a K_1,3 tree structure.
    c5 = 12
    
    # Print the coefficients in the specified order.
    print(f"c_1 = {c1}")
    print(f"c_2 = {c2}")
    print(f"c_3 = {c3}")
    print(f"c_4 = {c4}")
    print(f"c_5 = {c5}")

solve()