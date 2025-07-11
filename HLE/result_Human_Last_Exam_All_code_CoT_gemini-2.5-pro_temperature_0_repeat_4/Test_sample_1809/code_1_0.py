def solve_graph_orientation():
    """
    Determines the valid orientation number of the graph H.

    The solution is derived through logical deduction, as explained in the comments.
    """

    print("Step 1: Analyzing the constraints on vertex indegrees.")
    print("Let H be the graph. Let v1, v2, v3, v4 be the vertices of the initial K4.")
    print("Let U_i be the set of 30 vertices belonging to the 10 K3s attached to vi.")
    print("A valid orientation requires that for any adjacent vertices x and y, indegree(x) != indegree(y).\n")

    print("Step 2: Determining the possible indegrees for the attachment vertices (u-vertices).")
    print("Consider a single K3 attached to vi, with vertices {u1, u2, u3}.")
    print("These three vertices are mutually adjacent, and all are adjacent to vi.")
    print("Therefore, in a valid orientation, indeg(u1), indeg(u2), indeg(u3), and indeg(vi) must all be distinct.")
    print("The degree of any u-vertex is 3. By carefully orienting the edges within its K3 and the edge to vi, we can achieve distinct indegrees for {u1, u2, u3}.")
    print("The possible sets of indegrees for {u1, u2, u3} are subsets of {0, 1, 2, 3}.")
    print("For example:")
    print(" - If all 3 edges from the K3 point away from vi, the u-indegrees can be {1, 2, 3}.")
    print(" - If all 3 edges from the K3 point towards vi, the u-indegrees can be {0, 1, 2}.")
    print("In any case, for any vertex u in any attached K3, indeg(u) must be in the set {0, 1, 2, 3}.\n")

    print("Step 3: Determining the lower bound for the indegrees of the core vertices (v-vertices).")
    print("Each core vertex vi is adjacent to all 30 of its attached u-vertices in U_i.")
    print("Since indeg(u) is in {0, 1, 2, 3}, it must be that indeg(vi) is NOT in {0, 1, 2, 3}.")
    print("Therefore, for any valid orientation, indeg(vi) >= 4 for all i=1,2,3,4.")
    print("Furthermore, the four core vertices v1, v2, v3, v4 are mutually adjacent (they form a K4).")
    print("This means their indegrees must all be distinct.")
    print("The smallest possible set of four distinct integers, each at least 4, is {4, 5, 6, 7}.\n")

    print("Step 4: Establishing the lower bound for the valid orientation number.")
    print("From the above, any valid orientation must assign indegrees of at least 4, 5, 6, and 7 to the four core vertices.")
    print("This implies that the maximum indegree in any valid orientation must be at least 7.")
    print("So, the valid orientation number is >= 7.\n")

    print("Step 5: Constructing an orientation to show the number 7 is achievable.")
    print("We now show that a maximum indegree of 7 can be achieved. We need to find a valid orientation where max(indegree) = 7.")
    print("Let's assign the target indegrees {4, 5, 6, 7} to {v1, v2, v3, v4}.")
    print("The indegree of a core vertex is the sum of contributions from the K4 core and its attachments:")
    print("indeg(vi) = indeg_K4(vi) + indeg_attachments(vi)")
    print("First, orient the K4 core edges as a transitive tournament (vi -> vj for i < j).")
    print("This gives base indegrees: indeg_K4(v1)=0, indeg_K4(v2)=1, indeg_K4(v3)=2, indeg_K4(v4)=3.")
    print("Let C_i = indeg_attachments(vi). We can set up a system of equations:")
    
    # The final equations
    v1_k4, v2_k4, v3_k4, v4_k4 = 0, 1, 2, 3
    v1_target, v2_target, v3_target, v4_target = 4, 5, 6, 7
    c1 = v1_target - v1_k4
    c2 = v2_target - v2_k4
    c3 = v3_target - v3_k4
    c4 = v4_target - v4_k4

    print(f"  indeg(v1) = {v1_k4} + C1 = {v1_target}  => C1 = {c1}")
    print(f"  indeg(v2) = {v2_k4} + C2 = {v2_target}  => C2 = {c2}")
    print(f"  indeg(v3) = {v3_k4} + C3 = {v3_target}  => C3 = {c3}")
    print(f"  indeg(v4) = {v4_k4} + C4 = {v4_target}  => C4 = {c4}")

    print(f"\nThis requires an attachment contribution of C_i = {c1} for each core vertex.")
    print("This contribution is achievable. For each vi, we can orient its 10 K3s such that the total number of incoming edges is 4 (e.g., 4 K3s contribute 1 each, and 6 K3s contribute 0).")
    print("This construction yields a valid orientation where the v-indegrees are {4, 5, 6, 7} and all u-indegrees are in {0, 1, 2, 3}.")
    print("The maximum indegree in this orientation is 7.\n")

    final_answer = 7
    print("Conclusion: The smallest maximum indegree (the valid orientation number) is 7.")
    return final_answer

if __name__ == '__main__':
    solve_graph_orientation()
    # The final answer is printed as part of the explanation.
    # To conform to the required output format, we also print it here.
    print(f"\n<<<{7}>>>")