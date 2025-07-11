import math

def solve():
    """
    Solves the proportionality puzzle by calculating s1 and s2.
    """
    k = 100  # Committee size

    print("--- Step 1: Problem Setup and Definitions ---")
    print(f"Committee size k = {k}")
    print("Voter 1 must be unsatisfied, so candidates {a,b,c,x} are not in the committee W.")
    print("We want to find the smallest number of voters n such that a committee W exists that is 'fair' (PJR or EJR) but leaves voter 1 unsatisfied.")
    print("\nLet's analyze the conditions for PJR and EJR.")

    print("\n--- For Proportional Justified Representation (PJR) ---")
    print("A committee W fails PJR if there is a 1-cohesive group of unsatisfied voters S with size |S| >= n/k.")
    print("To satisfy PJR, for any such group S, we must have |S| < n/k, which means n > k * |S|.")
    print("To find the smallest n, we must choose W to minimize the size of the largest possible 1-cohesive group of unsatisfied voters, |S|_max.")
    print(f"Then, the minimum n (s1) must be floor(k * |S|_max) + 1.")

    print("\n--- For Extended Justified Representation (EJR) ---")
    print("A committee W fails EJR if for some l, there is an l-cohesive group of unsatisfied voters S with size |S| >= l*n/k.")
    print("To satisfy EJR, for ALL l and ALL such groups S, we must have |S| < l*n/k, which means n > k * |S| / l.")
    print("To find the smallest n, we must choose W to minimize the maximum value of (k * |S| / l) across all possible l and S.")
    print(f"Let this maximum be M. The minimum n (s2) must be floor(M) + 1.")
    
    print("\n--- Step 2: Analyzing Choices for the Committee W ---")
    print("Voter 1 ({a,b,c,x}) must be unsatisfied.")
    print("The satisfaction of voters 2-6 depends on whether y and z are in W.")
    print("Let's analyze four cases based on the inclusion of y and z in W.")

    # Data for the first 6 voters' cohesive groups
    # {group_name: (voter_indices, common_candidates)}
    groups_data = {
        'U1': ({1}, {'a', 'b', 'c', 'x'}),
        'U2_3': ({2, 3}, {'a', 'b', 'c', 'y'}),
        'U4_5_6': ({4, 5, 6}, {'a', 'b', 'c', 'z'}),
        'U1_to_3': ({1, 2, 3}, {'a', 'b', 'c'}),
        'U1_4_to_6': ({1, 4, 5, 6}, {'a', 'b', 'c'}),
        'U1_to_6': ({1, 2, 3, 4, 5, 6}, {'a', 'b', 'c'})
    }
    
    cases = [
        {"name": "Case 1: y not in W, z not in W", "unsatisfied": {1, 2, 3, 4, 5, 6}},
        {"name": "Case 2: y in W, z not in W", "unsatisfied": {1, 4, 5, 6}},
        {"name": "Case 3: z in W, y not in W", "unsatisfied": {1, 2, 3}},
        {"name": "Case 4: y in W, z in W", "unsatisfied": {1}}
    ]

    s1_min = float('inf')
    s2_min = float('inf')

    for case in cases:
        print(f"\n--- {case['name']} ---")
        unsatisfied_voters = case['unsatisfied']
        print(f"The set of unsatisfied voters (among the first 6) is: {sorted(list(unsatisfied_voters))}")

        max_s_size = 0
        max_ejr_bound = 0
        
        # We check all subgroups of the unsatisfied voters
        for group_name, (voter_indices, common_candidates) in groups_data.items():
            if voter_indices.issubset(unsatisfied_voters):
                s_size = len(voter_indices)
                cohesion_l = len(common_candidates)

                # PJR analysis
                max_s_size = max(max_s_size, s_size)

                # EJR analysis
                # For a group S, we need n > k * |S| / l for all l from 1 to L.
                # The tightest bound comes from l=1.
                ejr_bound = k * s_size / 1
                max_ejr_bound = max(max_ejr_bound, ejr_bound)

        # PJR calculation for this case
        pjr_bound = k * max_s_size
        n_pjr = math.floor(pjr_bound) + 1
        s1_min = min(s1_min, n_pjr)
        print(f"Largest unsatisfied cohesive group has size |S| = {max_s_size}.")
        print(f"PJR requires n > k * |S| = {k} * {max_s_size} = {pjr_bound}.")
        print(f"So, for this case, the smallest n for PJR is {n_pjr}.")
        
        # EJR calculation for this case
        n_ejr = math.floor(max_ejr_bound) + 1
        s2_min = min(s2_min, n_ejr)
        print(f"The maximum EJR bound is max(k * |S| / l). The tightest constraint for any group is when l=1.")
        print(f"Maximum bound max(k*|S|) = {max_ejr_bound}.")
        print(f"So, for this case, the smallest n for EJR is {n_ejr}.")


    print("\n--- Step 3: Final Results ---")
    print("To find the overall smallest n (s1 and s2), we take the minimum from all cases.")
    print("This corresponds to choosing the committee W that imposes the weakest constraint on n.")
    
    print(f"\nFor PJR, the minimum n (s1) is {s1_min}.")
    print(f"This is achieved in Case 4, where |S|_max = 1.")
    print(f"The condition is n > {k} * 1, which means n must be at least 101.")
    
    print(f"\nFor EJR, the minimum n (s2) is {s2_min}.")
    print(f"This is also achieved in Case 4, where the unsatisfied group is {{1}}.")
    print(f"The group {{1}} is l-cohesive for l=1,2,3,4. The conditions are n > {k}*1/l for each l.")
    print(f"The tightest condition is for l=1: n > {k}*1/1 = {k}.")
    print(f"n must be at least 101.")

    print("\nFinal Answer:")
    final_answer = (s1_min, s2_min)
    return final_answer
    
s1, s2 = solve()
print(f"s1 = {s1}")
print(f"s2 = {s2}")
print(f"<<<({s1}, {s2})>>>")
