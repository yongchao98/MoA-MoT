import math

def solve():
    """
    This script calculates the valid orientation number for the graph H.
    It systematically searches for the smallest possible maximum indegree (K)
    by checking if a valid orientation exists for K=1, 2, 3, ...
    """
    print("Determining the valid orientation number of graph H.")
    print("Let the four central vertices be v_0, v_1, v_2, v_3, with K_4-indegrees of 0, 1, 2, 3 respectively.")
    print("Let N_d be the number of K_3 groups with edges oriented towards v_d.")
    print("The total indegree of v_d is indeg(v_d) = d + 3 * N_d.\n")

    print("Constraints on N_d are:")
    print("1. 0 <= N_d <= 10")
    print("2. indeg(v_d) not in {0, 1, 2, 3} => N_0>=2, N_1>=1, N_2>=1, N_3>=1")
    print("3. indeg(v_0) != indeg(v_3) => N_0 != 1 + N_3\n")

    k = 1
    while True:
        # For a given max indegree k, find the allowed ranges for N_d
        n0_min, n1_min, n2_min, n3_min = 2, 1, 1, 1
        
        # Calculate max possible N_d values for a given k
        # indeg(v_d) <= k  => d + 3*N_d <= k => N_d <= (k-d)/3
        n0_max = min(math.floor(k / 3), 10)
        n1_max = min(math.floor((k - 1) / 3), 10)
        n2_max = min(math.floor((k - 2) / 3), 10)
        n3_max = min(math.floor((k - 3) / 3), 10)

        # Check if there is any valid range for each N_d
        if n0_min > n0_max or n1_min > n1_max or n2_min > n2_max or n3_min > n3_max:
            k += 1
            continue

        # Check if we can satisfy the N_0 != 1 + N_3 constraint
        found_solution = False
        solution_config = {}
        for n0_val in range(n0_min, n0_max + 1):
            for n3_val in range(n3_min, n3_max + 1):
                if n0_val != 1 + n3_val:
                    # A valid combination exists.
                    # We can pick the smallest valid N_1 and N_2.
                    solution_config = {'n0': n0_val, 'n1': n1_min, 'n2': n2_min, 'n3': n3_val}
                    found_solution = True
                    break
            if found_solution:
                break
        
        if found_solution:
            print(f"Found the smallest maximum indegree K = {k}.")
            print("A valid orientation can be constructed with this value.")
            
            n0 = solution_config['n0']
            n1 = solution_config['n1']
            n2 = solution_config['n2']
            n3 = solution_config['n3']

            print("\nExample configuration:")
            print(f"N_0 = {n0}, N_1 = {n1}, N_2 = {n2}, N_3 = {n3}")
            
            indeg_v0 = 0 + 3 * n0
            indeg_v1 = 1 + 3 * n1
            indeg_v2 = 2 + 3 * n2
            indeg_v3 = 3 + 3 * n3
            
            print("\nResulting indegrees of the central vertices:")
            print(f"indeg(v_0) = 0 + 3 * {n0} = {indeg_v0}")
            print(f"indeg(v_1) = 1 + 3 * {n1} = {indeg_v1}")
            print(f"indeg(v_2) = 2 + 3 * {n2} = {indeg_v2}")
            print(f"indeg(v_3) = 3 + 3 * {n3} = {indeg_v3}")
            
            max_indegree = max(indeg_v0, indeg_v1, indeg_v2, indeg_v3)
            print(f"\nThe maximum indegree in this orientation is {max_indegree}.")
            print("The indegrees of all other vertices are in {0, 1, 2, 3}.")
            print(f"Thus, the valid orientation number of H is {max_indegree}.")
            return max_indegree

        k += 1

final_answer = solve()
print(f"\n<<<9>>>")