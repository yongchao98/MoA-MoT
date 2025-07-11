import math

def solve_graph_orientation_problem():
    """
    Solves for the valid orientation number of the graph H.
    The logic is explained step-by-step.
    """
    
    print("### Step 1: Analyzing the Graph Structure and Definitions ###")
    print("The graph H is built from a central complete graph K_4 on vertices {v1, v2, v3, v4}.")
    print("To each vertex vi of the K_4, we attach 10 copies of K_3.")
    print("An attachment consists of connecting vi to all 3 vertices of a K_3 copy.")
    print("\nA 'valid orientation' of H assigns a direction to each edge such that for any adjacent vertices u and v, their indegrees d^-(u) and d^-(v) are different.")
    print("The 'valid orientation number' is the smallest possible maximum indegree over all valid orientations.")
    print("-" * 50)

    print("### Step 2: Deriving Constraints from the 'Valid Orientation' Rule ###")
    print("\n--- Part A: Constraints on vertices in a K_3 copy ---")
    print("Let's consider one K_3 copy attached to a central vertex, say v1. Let the K_3 vertices be u_a, u_b, u_c.")
    print("1. Inside the K_3, the vertices are mutually adjacent. Thus, d^-(u_a), d^-(u_b), and d^-(u_c) must all be different.")
    print("2. The orientation of edges within a K_3 must be a transitive tournament (e.g., a->b, b->c, a->c) to satisfy this. This results in internal indegrees of {0, 1, 2}.")
    print("3. Each of these vertices is also connected to v1. Let's analyze their final indegrees.")
    print("After re-labeling the K_3 vertices as u_0, u_1, u_2 according to their internal indegrees, their total indegrees are d^-(u_k) = k + delta_k, where delta_k is 0 or 1 depending on the orientation of the edge (v1, u_k).")
    print("To keep the indegrees {d^-(u_0), d^-(u_1), d^-(u_2)} distinct, only 4 of the 8 possible orientations of the edges to v1 are allowed. In all 4 allowed cases, the set of indegrees for {u_0, u_1, u_2} is a subset of {0, 1, 2, 3}.")

    print("\n--- Part B: Constraints on the central vertices v_i ---")
    print("A central vertex v_i is adjacent to all vertices in the K_3s attached to it.")
    print("Therefore, d^-(v_i) must be different from the indegree of any of these attached K_3 vertices.")
    print("Since the indegrees of the K_3 vertices are in the set {0, 1, 2, 3}, we must have d^-(v_i) not in {0, 1, 2, 3}.")
    min_indegree_v = 4
    print(f"This implies that for any central vertex v_i, its indegree must be at least {min_indegree_v}.")
    
    print("\n--- Part C: Constraints among the central vertices v_i ---")
    print("The vertices v1, v2, v3, v4 form a K_4, meaning they are all mutually adjacent.")
    print("Therefore, in any valid orientation, their indegrees d^-(v1), d^-(v2), d^-(v3), d^-(v4) must all be distinct from each other.")
    print("-" * 50)
    
    print("### Step 3: Calculating a Lower Bound for the Valid Orientation Number ###")
    print("From Step 2, we know the indegrees of the four central vertices must satisfy two conditions:")
    print("1. They must be four distinct integers.")
    print(f"2. Each must be greater than or equal to {min_indegree_v}.")
    print("\nTo find the minimum possible value for the maximum of these four indegrees, we should choose the smallest four distinct integers that satisfy condition 2.")
    
    required_indegrees = []
    current_val = min_indegree_v
    for i in range(4):
        required_indegrees.append(current_val)
        current_val += 1
        
    lower_bound = max(required_indegrees)
    print(f"The smallest possible set of indegrees for {{v1, v2, v3, v4}} is {{{', '.join(map(str, required_indegrees))}}}.")
    print(f"This means the maximum indegree in any valid orientation must be at least {lower_bound}.")
    print(f"So, the valid orientation number is >= {lower_bound}.")
    print("-" * 50)

    print("### Step 4: Constructing an Orientation to Prove the Lower Bound is Achievable ###")
    print("We will now show that a maximum indegree of 7 is achievable.")
    print("The indegree of a vertex v_i is d^-(v_i) = d^-(v_i)_from_K4 + d^-(v_i)_from_K3s.")
    
    print("\nFirst, orient the central K_4. For the 4 vertices to have distinct indegrees from their 6 internal edges, their indegrees must be {0, 1, 2, 3}.")
    k4_indegrees = {'v1': 0, 'v2': 1, 'v3': 2, 'v4': 3}
    
    print(f"Let's assign these indegrees: d^-(v1)_K4 = {k4_indegrees['v1']}, d^-(v2)_K4 = {k4_indegrees['v2']}, d^-(v3)_K4 = {k4_indegrees['v3']}, d^-(v4)_K4 = {k4_indegrees['v4']}.")
    
    print("\nNext, orient the edges connecting the v_i to their K_3s. Let S_i be the sum of indegree contributions from the 10 K_3s attached to v_i.")
    print("We aim to make the final indegrees for the v_i vertices be {4, 5, 6, 7}. We solve for each S_i:")
    
    s_values = {}
    print("\n--- The Final Indegree Equations ---")
    for i in range(1, 5):
        vi_name = f'v{i}'
        s_values[vi_name] = required_indegrees[i-1] - k4_indegrees[vi_name]
        print(f"d^-({vi_name}) = d^-({vi_name})_K4 + S{i}   =>   {required_indegrees[i-1]} = {k4_indegrees[vi_name]} + S{i}")
        # Note: I swapped my target assignment to make all Si the same. Let's fix that.
        # d^-(v1) = 0 + S1, d^-(v2) = 1 + S2, etc. Let's assign targets to them.
        # Target d^-(v1)=4, d^-(v2)=5, d^-(v3)=6, d^-(v4)=7.
        # S1=4, S2=4, S3=4, S4=4.
    
    common_s_value = 4
    print(f"\nWe can achieve the target indegrees {required_indegrees} by setting each S_i = {common_s_value}.")
    
    final_v_indegrees = {}
    for i in range(1, 5):
        vi_name = f'v{i}'
        final_v_indegrees[vi_name] = k4_indegrees[vi_name] + common_s_value
        print(f"d^-({vi_name}) = {k4_indegrees[vi_name]} + {common_s_value} = {final_v_indegrees[vi_name]}")
        
    print(f"\nAchieving S_i = {common_s_value} for each v_i is possible. For example, for each v_i, orient the edges to its 10 K_3s such that two K_3 groups contribute 2 each to the indegree, and the other eight contribute 0.")
    print("This construction results in a valid orientation with v_i indegrees {4, 5, 6, 7}.")
    
    max_u_indegree = 3
    achieved_max_indegree = max(final_v_indegrees.values())
    overall_max_indegree = max(achieved_max_indegree, max_u_indegree)
    
    print(f"\nThe maximum indegree among the v_i is {achieved_max_indegree}.")
    print(f"The maximum possible indegree for any u-vertex is {max_u_indegree}.")
    print(f"The overall maximum indegree in this orientation is max({achieved_max_indegree}, {max_u_indegree}) = {overall_max_indegree}.")
    print(f"Thus, the valid orientation number is <= {overall_max_indegree}.")
    print("-" * 50)
    
    print("### Step 5: Conclusion ###")
    print(f"We have shown that the valid orientation number is >= {lower_bound}.")
    print(f"We have also shown that the valid orientation number is <= {overall_max_indegree}.")
    print(f"Since {lower_bound} == {overall_max_indegree}, we can conclude the final answer.")

    final_answer = overall_max_indegree
    print(f"\nThe valid orientation number of graph H is {final_answer}.")
    return final_answer

solve_graph_orientation_problem()