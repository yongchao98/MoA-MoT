def solve_graph_orientation():
    """
    Calculates the valid orientation number of the graph H.

    The code follows a logical deduction to find the smallest possible maximum indegree
    in a valid orientation of the graph H.
    """

    # Step 1: Define the problem and constraints.
    # The graph H has 4 central vertices (v_1, v_2, v_3, v_4) and for each v_i,
    # 10 K_3s are attached. A valid orientation requires adjacent vertices to have
    # different indegrees. We want to find the smallest possible maximum indegree.

    # Step 2: Analyze the indegrees of central vertices (d_i).
    # The four central vertices v_i are all adjacent, so their indegrees d_1, d_2, d_3, d_4
    # must be distinct.
    # A central vertex v_i is also adjacent to 30 peripheral vertices. Its indegree d_i
    # must be different from the indegrees of all its neighbors.

    # Step 3: Develop an optimal orientation strategy.
    # The indegree of a central vertex v_i is d_i = i_i + c_i, where:
    # - i_i is its indegree from the central K_4.
    # - c_i is its indegree from its 30 peripheral neighbors.

    # To minimize the maximum d_i, we make optimal choices for i_i and c_i.

    # Choice for central K_4 indegrees (i_i):
    # We orient the K_4 edges acyclically (e.g., v_i -> v_j for i < j).
    # This gives the smallest possible distinct non-negative indegrees.
    i = [0, 1, 2, 3]
    print(f"To minimize indegrees, we orient the central K_4 to give indegrees: {i[0]}, {i[1]}, {i[2]}, {i[3]}.")

    # Choice for peripheral connection indegrees (c_i):
    # For each v_i, we must orient its 30 connecting edges. A careful analysis shows that
    # to ensure the indegrees of peripheral neighbors are distinct among themselves and
    # that their set of values doesn't conflict with d_i, the number of edges
    # pointing towards v_i (c_i) must be 0, 10, 20, or 30.

    # To minimize d_i = i_i + c_i, we must choose the smallest possible valid c_i.
    # A choice of c_i = 0 is only valid if i_i = 0.
    # For i_i > 0, the smallest valid choice for c_i is 10.

    # Step 4: Calculate the indegrees in the optimal orientation.

    # For v_1, with i_1 = 0:
    i_1 = 0
    c_1 = 0
    d_1 = i_1 + c_1
    print(f"\nFor central vertex v_1, K_4 indegree is {i_1}. We can choose peripheral indegree contribution c_1 = {c_1}.")
    print(f"Indegree of v_1 is d_1 = {i_1} + {c_1} = {d_1}.")

    # For v_2, with i_2 = 1:
    i_2 = 1
    c_2 = 10
    d_2 = i_2 + c_2
    print(f"\nFor central vertex v_2, K_4 indegree is {i_2}. The minimum valid peripheral contribution is c_2 = {c_2}.")
    print(f"Indegree of v_2 is d_2 = {i_2} + {c_2} = {d_2}.")

    # For v_3, with i_3 = 2:
    i_3 = 2
    c_3 = 10
    d_3 = i_3 + c_3
    print(f"\nFor central vertex v_3, K_4 indegree is {i_3}. The minimum valid peripheral contribution is c_3 = {c_3}.")
    print(f"Indegree of v_3 is d_3 = {i_3} + {c_3} = {d_3}.")

    # For v_4, with i_4 = 3:
    i_4 = 3
    c_4 = 10
    d_4 = i_4 + c_4
    print(f"\nFor central vertex v_4, K_4 indegree is {i_4}. The minimum valid peripheral contribution is c_4 = {c_4}.")
    print(f"Indegree of v_4 is d_4 = {i_4} + {c_4} = {d_4}.")

    # Step 5: Determine the maximum indegree.
    central_indegrees = [d_1, d_2, d_3, d_4]
    max_central_indegree = max(central_indegrees)

    # The indegrees of peripheral vertices in this construction are all in {0, 1, 2, 3}.
    max_peripheral_indegree = 3
    print(f"\nThe maximum indegree for any peripheral vertex in this orientation is {max_peripheral_indegree}.")

    valid_orientation_number = max(max_central_indegree, max_peripheral_indegree)
    print(f"\nThe set of central indegrees is {{{d_1}, {d_2}, {d_3}, {d_4}}}.")
    print(f"The maximum indegree in the graph is max({max_central_indegree}, {max_peripheral_indegree}) = {valid_orientation_number}.")
    print(f"\nThe valid orientation number of H is {valid_orientation_number}.")


solve_graph_orientation()
<<<13>>>