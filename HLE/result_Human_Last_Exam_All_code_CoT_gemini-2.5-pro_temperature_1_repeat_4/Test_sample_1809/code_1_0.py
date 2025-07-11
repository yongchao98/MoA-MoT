def solve_graph_orientation_problem():
    """
    This script calculates the valid orientation number of the graph H.
    The method is based on a step-by-step logical deduction, as explained above.
    """
    print("Method to find the valid orientation number of graph H:")
    print("1. The graph H has 4 central vertices (v1, v2, v3, v4) forming a K4.")
    print("2. Each central vertex vi has 10 K3s attached. For any K3, its 3 vertices and vi form a K4 subgraph.")
    print("3. In a valid orientation, the indegrees of all 4 vertices in such a K4 subgraph must be distinct.")
    print("4. The indegrees of the K3 vertices can only be 0, 1, 2, or 3.")
    print("5. This implies that for each central vertex vi, the indegrees of all its 30 peripheral neighbors must come from the same set S_i of 3 values from {0,1,2,3}.")
    print("6. The indegree of a central vertex, d(vi), is the sum of its indegree from the central K4 (ci) and from its peripheral neighbors (dpi).")
    print("7. The values for {c1, c2, c3, c4} must be a permutation of {0, 1, 2, 3}.")
    print("8. The values for dpi are determined by the choice of S_i and must be from {0, 10, 20, 30}.")
    print("9. A key constraint arises: if dpi = 0, then ci must be 0. This means at most one vertex can have dpi=0.")
    print("\nTo find the valid orientation number, we must find an assignment that minimizes the maximum indegree.")
    print("Let's construct an optimal orientation:")
    print("- Assign dpi = 10 for all four central vertices (v1, v2, v3, v4).")
    print("- This corresponds to choosing the peripheral indegree set S_i = {0, 2, 3} for all i.")
    print("- Assign the central K4 indegrees ci as a permutation of {0, 1, 2, 3}.")
    print("- The total indegrees of the central vertices d(vi) = ci + dpi are:")
    
    c_values = [0, 1, 2, 3]
    dp_value = 10
    
    d_v1 = c_values[0] + dp_value
    d_v2 = c_values[1] + dp_value
    d_v3 = c_values[2] + dp_value
    d_v4 = c_values[3] + dp_value
    
    print(f"  d(v1) = {c_values[0]} + {dp_value} = {d_v1}")
    print(f"  d(v2) = {c_values[1]} + {dp_value} = {d_v2}")
    print(f"  d(v3) = {c_values[2]} + {dp_value} = {d_v3}")
    print(f"  d(v4) = {c_values[3]} + {dp_value} = {d_v4}")
    
    max_indegree = max(d_v1, d_v2, d_v3, d_v4)
    
    print(f"\nThe indegrees of peripheral vertices are in {{0, 2, 3}}, so their maximum is 3.")
    print(f"The maximum indegree in this orientation is max({d_v1}, {d_v2}, {d_v3}, {d_v4}, 3) = {max_indegree}.")
    print("This construction is valid and achieves a maximum indegree of 13.")
    print("Any other choice of dpi values results in a maximum indegree of 13 or greater.")
    print("\nTherefore, the valid orientation number of H is 13.")

solve_graph_orientation_problem()