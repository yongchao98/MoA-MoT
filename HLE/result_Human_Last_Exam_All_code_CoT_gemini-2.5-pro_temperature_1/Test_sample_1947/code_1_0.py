def solve_coefficients():
    """
    This function determines the coefficients for the expression of the number of
    closed tree-like walks of length 6 in a simple graph X.

    The expression is of the form:
    c1*e + c2*k + c3*p + c4*sum(deg(v) choose 2) + c5*sum(deg(v) choose 3)

    Based on combinatorial counting for each possible underlying tree structure of the walk,
    the coefficients are determined as follows:
    - c1 corresponds to walks on a P_2 graph (a single edge).
    - c2 corresponds to walks on a K_3 graph (a triangle), which is not a tree.
    - c3 corresponds to walks on a P_4 graph (a path of length 3).
    - c4 corresponds to walks on a P_3 graph (a path of length 2).
    - c5 corresponds to walks on a K_1,3 graph (a claw).
    """

    # c1: Underlying graph is a P_2 (1 edge), traversed 6 times.
    # Walks: u-v-u-v-u-v-u and v-u-v-u-v-u-v. 2 walks per edge.
    c1 = 2

    # c2: Underlying graph cannot be a K_3 (triangle) as it's not a tree.
    c2 = 0

    # c3: Underlying graph is a P_4 (3 edges), each traversed twice.
    # 6 walks per P_4 subgraph.
    c3 = 6

    # c4: Underlying graph is a P_3 (2 edges), one traversed 4 times, one twice.
    # 12 walks per P_3 subgraph.
    c4 = 12

    # c5: Underlying graph is a K_1,3 (3 edges), each traversed twice.
    # 12 walks per K_1,3 subgraph.
    c5 = 12

    coefficients = [c1, c2, c3, c4, c5]
    
    # "Remember in the final code you still need to output each number in the final equation!"
    # The numbers in the final equation are the coefficients c_1, ..., c_5.
    for c in coefficients:
        print(c)

solve_coefficients()