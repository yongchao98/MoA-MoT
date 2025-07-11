def solve_graph_orientation_problem():
    """
    This script derives the valid orientation number for the graph H.
    The derivation follows a logical, step-by-step process.
    """
    
    print("### Step 1 & 2: Analyze Graph and Constraints ###")
    print("The graph H has a central complete graph K_4 (vertices v1, v2, v3, v4).")
    print("Each vertex v_i is connected to 10 copies of K_3.")
    print("A valid orientation requires adjacent vertices to have different indegrees.")
    print("-" * 30)

    print("### Step 3: Devise an Orientation Strategy ###")
    print("To ensure vertices within complete subgraphs have distinct indegrees, we orient them transitively.")
    print("For K_4 (v1, v2, v3, v4), we orient edges v_i -> v_j if i < j.")
    print("This gives base indegrees from within the K_4:")
    k4_base_indegrees = {'v1': 0, 'v2': 1, 'v3': 2, 'v4': 3}
    for v, d in k4_base_indegrees.items():
        print(f"  - Base indegree of {v}: {d}")
    print("\nFor each K_3, we orient it transitively, giving its vertices base indegrees of 0, 1, and 2.")
    print("-" * 30)

    print("### Step 4: Analyze Indegrees ###")
    print("The indegree of a K_3 vertex is its base indegree (0, 1, or 2) plus a contribution from its edge to the central v_i (0 or 1).")
    print("Thus, the indegree of any K_3 vertex must be in the set {0, 1, 2, 3}.")
    
    print("\nThe indegree of a central vertex v_i is its base indegree plus the contribution from its 30 connecting edges to the triangles.")
    print("Let C_i be this contribution for v_i. d(v_i) = base_indegree(v_i) + C_i.")
    print("\nConstraints on d(v_i):")
    print("1. The four values d(v1), d(v2), d(v3), d(v4) must be distinct because the v_i are all mutually adjacent.")
    print("2. d(v_i) must be different from the indegree of any adjacent K_3 vertex. So, d(v_i) cannot be 0, 1, 2, or 3.")
    print("-" * 30)

    print("### Step 5: Calculate the Lower Bound ###")
    print("We need to find four distinct integer indegrees for the v_i vertices, each of which must be 4 or greater.")
    print("The smallest possible set of such integers is {4, 5, 6, 7}.")
    smallest_set_of_indegrees = [4, 5, 6, 7]
    print(f"Smallest possible indegrees for v1, v2, v3, v4: {smallest_set_of_indegrees}")
    
    lower_bound = max(smallest_set_of_indegrees)
    print(f"\nThe maximum indegree in any valid orientation must be at least the maximum of this set.")
    print(f"So, the valid orientation number must be >= {lower_bound}.")
    print("-" * 30)

    print("### Step 6: Show Achievability ###")
    print(f"We now show that a maximum indegree of {lower_bound} is achievable.")
    print("We assign the target indegrees {4, 5, 6, 7} to v1, v2, v3, v4 respectively.")
    
    target_d_v1 = 4
    C1 = target_d_v1 - k4_base_indegrees['v1']
    print(f"d(v1) = base(v1) + C1 => {target_d_v1} = {k4_base_indegrees['v1']} + C1 => C1 = {C1}")

    target_d_v2 = 5
    C2 = target_d_v2 - k4_base_indegrees['v2']
    print(f"d(v2) = base(v2) + C2 => {target_d_v2} = {k4_base_indegrees['v2']} + C2 => C2 = {C2}")

    target_d_v3 = 6
    C3 = target_d_v3 - k4_base_indegrees['v3']
    print(f"d(v3) = base(v3) + C3 => {target_d_v3} = {k4_base_indegrees['v3']} + C3 => C3 = {C3}")

    target_d_v4 = 7
    C4 = target_d_v4 - k4_base_indegrees['v4']
    print(f"d(v4) = base(v4) + C4 => {target_d_v4} = {k4_base_indegrees['v4']} + C4 => C4 = {C4}")

    print("\nFor each v_i, a contribution of C_i=4 is needed from its 30 connecting edges.")
    print("This is feasible: for each v_i, orient 4 of these 30 edges towards v_i and 26 away.")
    print(f"This construction results in a valid orientation with maximum indegree {lower_bound}.")
    print("-" * 30)

    print("### Step 7: Final Conclusion ###")
    valid_orientation_number = lower_bound
    print(f"The valid orientation number is >= {lower_bound} and <= {lower_bound}.")
    print(f"Therefore, the valid orientation number of H is {valid_orientation_number}.")
    
    return valid_orientation_number

if __name__ == '__main__':
    solve_graph_orientation_problem()
