def solve_graph_counting():
    """
    Calculates the numbers n_j of isomorphism classes of vertex-transitive graphs
    with 8 vertices and vertex degree precisely j for j=0, ..., 7.
    The script follows a step-by-step logical deduction based on graph theory principles.
    """
    print("Calculating the number of 8-vertex vertex-transitive graphs by degree j...")
    
    n_counts = [0] * 8

    # Step 1: Degree j=0
    # The only 0-regular graph on 8 vertices is the empty graph. It is unique
    # and vertex-transitive.
    n_0 = 1
    n_counts[0] = n_0
    print(f"For degree j=0, the graph is the empty graph. It is vertex-transitive. So, n_0 = {n_0}")

    # Step 2: Degree j=1
    # The only 1-regular graph on 8 vertices is the perfect matching (4 disjoint edges).
    # It is unique up to isomorphism and is vertex-transitive.
    n_1 = 1
    n_counts[1] = n_1
    print(f"For degree j=1, the graph is a perfect matching (4K2). It is vertex-transitive. So, n_1 = {n_1}")

    # Step 3: Degree j=2
    # 2-regular graphs are disjoint unions of cycles. On 8 vertices, the partitions
    # of 8 into parts of size >= 3 are [8], [5, 3], and [4, 4].
    # - C_8: Vertex-transitive.
    # - C_5 + C_3: Not vertex-transitive (vertices in C_5 cannot be mapped to vertices in C_3).
    # - C_4 + C_4: Vertex-transitive (automorphisms can swap the two C_4 components).
    # This gives 2 graphs.
    n_2 = 2
    n_counts[2] = n_2
    print(f"For degree j=2, the graphs are C_8 and 2C_4. Both are vertex-transitive. So, n_2 = {n_2}")

    # Step 4: Degree j=3
    # There are 6 non-isomorphic 3-regular (cubic) graphs on 8 vertices.
    # It is a known result that all 6 of these graphs are vertex-transitive.
    n_3 = 6
    n_counts[3] = n_3
    print(f"For degree j=3, all 6 cubic graphs on 8 vertices are vertex-transitive. So, n_3 = {n_3}")

    # Step 5: Use symmetry for degrees j=4, 5, 6, 7
    # The number of vertex-transitive graphs of degree j is the same as for degree 7-j.
    n_4 = n_3
    n_counts[4] = n_4
    print(f"For degree j=4, by complement symmetry, n_4 = n_3 = {n_4}")

    n_5 = n_2
    n_counts[5] = n_5
    print(f"For degree j=5, by complement symmetry, n_5 = n_2 = {n_5}")

    n_6 = n_1
    n_counts[6] = n_6
    print(f"For degree j=6, by complement symmetry, n_6 = n_1 = {n_6}")

    n_7 = n_0
    n_counts[7] = n_7
    print(f"For degree j=7, by complement symmetry, n_7 = n_0 = {n_7}")

    # Step 6: Print the final result in the specified format
    result_list_str = f"[{n_counts[0]}, {n_counts[1]}, {n_counts[2]}, {n_counts[3]}, {n_counts[4]}, {n_counts[5]}, {n_counts[6]}, {n_counts[7]}]"
    
    print("\n" + "="*40)
    print("The final list of counts [n_0, n_1, ..., n_7] is:")
    print(result_list_str)
    print("="*40)

solve_graph_counting()