def solve():
    """
    This function determines the coefficients for the expression of the number of
    closed tree-like walks of length 6 in a simple graph X.

    The expression is of the form:
    c_1 * e + c_2 * k + c_3 * p + c_4 * sum(deg(v) choose 2) + c_5 * sum(deg(v) choose 3)

    where:
    - e is the number of edges.
    - k is the number of K_3 subgraphs (triangles).
    - p is the number of P_4 subgraphs (paths of length 3).
    - deg(v) is the degree of vertex v.
    """

    # c_1: Contribution from P_2 trees (single edges)
    # An edge (u,v) is traversed 6 times. Walks: u->v->u->... or v->u->v->...
    # 2 walks per edge.
    c1 = 2

    # c_2: Contribution from K_3 subgraphs
    # A tree-like walk's edges must form a tree. A K_3 is a cycle.
    # The count of walks on a specific tree structure is independent of whether
    # that tree is part of a larger cycle (like a triangle).
    # Thus, the k term is not needed.
    c2 = 0

    # c_3: Contribution from P_4 trees (paths of 3 edges)
    # Each of the 3 edges is traversed twice.
    # For a path v0-v1-v2-v3, we count walks starting at each vertex.
    # Start v0: 1, Start v1: 2, Start v2: 2, Start v3: 1. Total = 6 walks per P_4.
    c3 = 6

    # c_4: Contribution from P_3 trees (paths of 2 edges)
    # The 2 edges are traversed (4, 2) or (2, 4) times.
    # For a path u-v-w, total walks are 12.
    c4 = 12

    # c_5: Contribution from K_1,3 trees (star graphs with 3 edges)
    # Each of the 3 edges is traversed twice.
    # For a star with center v0 and leaves v1,v2,v3:
    # Start v0: 3! = 6. Start v1,v2,v3: 2 each. Total = 6+2+2+2 = 12 walks per K_1,3.
    c5 = 12

    # The final coefficients
    coefficients = [c1, c2, c3, c4, c5]
    
    # Printing each coefficient for the final equation
    print(f"c_1 = {coefficients[0]}")
    print(f"c_2 = {coefficients[1]}")
    print(f"c_3 = {coefficients[2]}")
    print(f"c_4 = {coefficients[3]}")
    print(f"c_5 = {coefficients[4]}")
    
    # The final answer in the required format
    answer = ",".join(map(str, coefficients))
    print(f"\nThe coefficients in order are: {answer}")

solve()