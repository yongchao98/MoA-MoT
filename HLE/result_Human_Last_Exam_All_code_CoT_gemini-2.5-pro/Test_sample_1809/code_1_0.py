import sys

def solve():
    """
    This function explains the step-by-step derivation of the valid orientation number for graph H.
    """
    
    # Step 1: Analyze the graph structure H
    print("Step 1: Analyzing the structure of graph H.")
    print("The graph H is constructed from a complete graph K_4 (on vertices v_1, v_2, v_3, v_4).")
    print("For each vertex v_i in K_4, 10 disjoint copies of K_3 are attached.")
    print("The connection is made by adding an edge from v_i to all 3 vertices of each of its 10 associated K_3s.\n")

    # Step 2: Calculate vertex degrees
    print("Step 2: Calculating the degrees of vertices in H.")
    # A core vertex v_i is connected to 3 other core vertices and 10 K_3s (each with 3 vertices).
    deg_vi = 3 + 10 * 3
    print(f"A core vertex v_i is connected to 3 other core vertices and 10 * 3 = 30 triangle vertices.")
    print(f"Therefore, the degree is deg(v_i) = 3 + 30 = {deg_vi}.\n")

    # A triangle vertex u is connected to 2 other vertices in its own K_3 and 1 core vertex.
    deg_u = 2 + 1
    print(f"A triangle vertex u is connected to 2 other vertices in its K_3 and 1 core vertex v_i.")
    print(f"Therefore, the degree is deg(u) = 2 + 1 = {deg_u}.\n")

    # Step 3: Lower bound for the valid orientation number (k)
    print("Step 3: Finding a lower bound for the valid orientation number, k.")
    print("A valid orientation requires that any two adjacent vertices have different indegrees.")
    
    print("\nConstraint 1: The four core vertices v_1, v_2, v_3, v_4 are all mutually adjacent.")
    print("This means their indegrees, let's call them d_1, d_2, d_3, d_4, must all be distinct from each other.\n")

    print("Constraint 2: A triangle vertex u has degree 3, so its indegree can only be 0, 1, 2, or 3.")
    print("The three vertices of any single K_3 are mutually adjacent, so their indegrees must be distinct.")
    print("By orienting the K_3's internal edges transitively (e.g., u1->u2, u1->u3, u2->u3), the indegrees from within the K_3 are 0, 1, and 2.")
    print("The edge to the core vertex v_i can add at most 1 to each of these indegrees.")
    print("This ensures that the final indegree of any triangle vertex u is always in the set {0, 1, 2, 3}.\n")

    print("Constraint 3: A triangle vertex u (attached to v_i) is adjacent to v_i.")
    print("Therefore, indeg(v_i) must be different from indeg(u).")
    print("Since indeg(u) must be in {0, 1, 2, 3}, the indegree of a core vertex, d_i = indeg(v_i), cannot be 0, 1, 2, or 3.")
    print(f"Thus, the indegree of any core vertex v_i must be at least 4.\n")

    print("Combining these constraints:")
    print("We have four distinct core indegrees d_1, d_2, d_3, d_4, and each must be at least 4.")
    print("The smallest possible set of such distinct integers is {4, 5, 6, 7}.")
    d_min_max = 7
    print(f"So, in any valid orientation, the maximum indegree among the core vertices must be at least {d_min_max}.")
    print(f"The maximum indegree of the whole graph, k, must be at least this value.")
    print(f"This establishes the lower bound: k >= {d_min_max}.\n")

    # Step 4: Upper bound for the valid orientation number (k)
    print("Step 4: Finding an upper bound for k by constructing a specific valid orientation.")
    print("We will construct a valid orientation with a maximum indegree of 7.")
    
    print("\nConstruction part 1: Orient the base K_4 transitively (e.g., v_i -> v_j if i < j).")
    indeg_k4 = [0, 1, 2, 3]
    print(f"This gives indegrees from K_4 edges as {indeg_k4[0]}, {indeg_k4[1]}, {indeg_k4[2]}, and {indeg_k4[3]} for v_1, v_2, v_3, v_4 respectively.\n")

    print("Construction part 2: Orient the edges between core vertices and their attached triangles.")
    print("Our goal is to make the final indegrees of v_1, v_2, v_3, v_4 equal to 4, 5, 6, and 7.")
    target_indeg_v = [4, 5, 6, 7]
    needed_in_edges = 4
    print(f"To get indeg(v_1) = {target_indeg_v[0]}, we need {target_indeg_v[0]} - {indeg_k4[0]} = {needed_in_edges} incoming edges from its 10 K_3s.")
    print(f"To get indeg(v_2) = {target_indeg_v[1]}, we need {target_indeg_v[1]} - {indeg_k4[1]} = {needed_in_edges} incoming edges.")
    print(f"To get indeg(v_3) = {target_indeg_v[2]}, we need {target_indeg_v[2]} - {indeg_k4[2]} = {needed_in_edges} incoming edges.")
    print(f"To get indeg(v_4) = {target_indeg_v[3]}, we need {target_indeg_v[3]} - {indeg_k4[3]} = {needed_in_edges} incoming edges.\n")

    print(f"Construction part 3: Show that achieving {needed_in_edges} incoming edges from 10 K_3s is possible.")
    print("For each K_3 attached to v_i, we can orient its edges to contribute 0, 1, 2, or 3 incoming edges to v_i, while keeping the K_3 vertices' indegrees valid.")
    n_k3_3 = 1
    n_k3_1 = 1
    n_k3_0 = 8
    total_in_edges = n_k3_3 * 3 + n_k3_1 * 1 + n_k3_0 * 0
    print(f"For each v_i, we can orient its 10 K_3s as follows:")
    print(f"- {n_k3_3} K_3 contributes 3 edges to v_i (indegrees of its vertices become {{0,1,2}}).")
    print(f"- {n_k3_1} K_3 contributes 1 edge to v_i (indegrees of its vertices become {{0,2,3}}).")
    print(f"- {n_k3_0} K_3s contribute 0 edges to v_i (indegrees of their vertices become {{1,2,3}}).")
    print(f"The total number of incoming edges is {n_k3_3} * 3 + {n_k3_1} * 1 + {n_k3_0} * 0 = {total_in_edges}, as required.\n")

    print("Construction part 4: Verify the constructed orientation.")
    print(f"- Indegrees of core vertices v_i are {target_indeg_v[0]}, {target_indeg_v[1]}, {target_indeg_v[2]}, {target_indeg_v[3]}. They are distinct, as required.")
    print("- Indegrees of all triangle vertices u are in the set {0, 1, 2, 3}.")
    print("- All adjacency conditions are met:")
    print(f"  - indeg(v_i) != indeg(v_j) (e.g., {target_indeg_v[0]} != {target_indeg_v[1]}).")
    print(f"  - indeg(v_i) != indeg(u) (e.g., {target_indeg_v[0]} is not in {{0,1,2,3}}).")
    print("  - indeg(u) != indeg(w) for u, w in the same K_3 (ensured by construction).\n")

    max_indeg = target_indeg_v[3]
    print(f"The maximum indegree in this orientation is max({target_indeg_v[3]}, 3) = {max_indeg}.")
    print(f"This establishes the upper bound: k <= {max_indeg}.\n")

    # Step 5: Conclusion
    print("Step 5: Conclusion.")
    print(f"From the lower bound (k >= {d_min_max}) and the upper bound (k <= {max_indeg}), we conclude the final answer.")
    final_answer = 7
    print(f"\nThe valid orientation number of H is {final_answer}.")
    
    return final_answer

if __name__ == "__main__":
    solve()
    # The final answer is wrapped according to the required format.
    # The code execution above prints the reasoning.
    # sys.stdout.write("<<<7>>>") # This would be for programmatic extraction.
                                 # We will just print the final result clearly.