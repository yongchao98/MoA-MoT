def solve_ejr_committee_problem():
    """
    Solves the EJR committee problem by defining the constraints and then
    constructing committees to find the minimum and maximum representation for Voter 1.
    """
    # --- 1. Define Problem Parameters ---
    n_voters = 10
    k_committee_size = 5
    quota = n_voters / k_committee_size
    
    print("--- Problem Setup ---")
    print(f"Number of voters (n): {n_voters}")
    print(f"Committee size (k): {k_committee_size}")
    print(f"Quota (n/k): {quota}")
    print("-" * 20 + "\n")

    # --- 2. Define Ballots and Key Sets ---
    voter_ballots = {
        'V1': {'x1','x2','x3','y1','z3'},
        'V2': {'x1','x2','x3','y2'},
        'V3': {'x1','x2','x3','y3'},
        'V4': {'y4','y5','y6','z1'},
        'V5': {'y4','y5','y6','z1'},
        'V6': {'y4','y5','y6','z2'},
        'V7': {'y4','y5','y6','z2'},
        'V8': {'x1','x2','x3','z1'},
        'V9': {'x1','x2','x3','z1'},
        'V10': {'x1','x2','x3','z1'}
    }
    C_V1 = voter_ballots['V1']

    # --- 3. Analyze Cohesive Groups and EJR Constraints ---
    # Group X: {V1, V2, V3, V8, V9, V10} -> 3-cohesive
    group_x_voters = ['V1', 'V2', 'V3', 'V8', 'V9', 'V10']
    U_X = set.union(*(voter_ballots[v] for v in group_x_voters))
    l_X = 3
    
    # Group Y: {V4, V5, V6, V7} -> 2-cohesive
    group_y_voters = ['V4', 'V5', 'V6', 'V7']
    U_Y = set.union(*(voter_ballots[v] for v in group_y_voters))
    l_Y = 2

    print("--- EJR Constraints ---")
    print(f"Group X ({len(group_x_voters)} voters) is {l_X}-cohesive.")
    print(f"Constraint 1: Committee must have at least {l_X} members from {sorted(list(U_X))}")
    print(f"Group Y ({len(group_y_voters)} voters) is {l_Y}-cohesive.")
    print(f"Constraint 2: Committee must have at least {l_Y} members from {sorted(list(U_Y))}")
    print("-" * 20 + "\n")

    # --- 4. Finding the Maximum ---
    print("--- Finding the Maximum Representation for Voter 1 ---")
    # Proof: |W| >= |W ∩ C_V1| + |W ∩ U_Y| => 5 >= m + 2 => m <= 3
    max_theoretical = k_committee_size - l_Y
    print(f"The number of V1-approved candidates is at most {k_committee_size} - {l_Y} = {max_theoretical}.")

    # Construct a committee W_max to demonstrate max is 3
    W_max = {'x1', 'x2', 'x3', 'y4', 'y5'}
    print(f"Let's test the committee W_max = {sorted(list(W_max))}")
    
    # Check constraints
    check1_max = len(W_max.intersection(U_X))
    check2_max = len(W_max.intersection(U_Y))
    print(f"Check 1: |W_max ∩ U_X| = {check1_max} (must be >= {l_X}) -> {'OK' if check1_max >= l_X else 'FAIL'}")
    print(f"Check 2: |W_max ∩ U_Y| = {check2_max} (must be >= {l_Y}) -> {'OK' if check2_max >= l_Y else 'FAIL'}")
    
    # Calculate representation for Voter 1
    max_approved = len(W_max.intersection(C_V1))
    print(f"Number of candidates approved by Voter 1: |W_max ∩ C_V1| = {len(W_max.intersection(C_V1))}")
    print(f"Final Maximum: {max_approved}\n")


    # --- 5. Finding the Minimum ---
    print("--- Finding the Minimum Representation for Voter 1 ---")
    
    # Construct a committee W_min to aim for min=0
    W_min = {'y2', 'y3', 'z1', 'y4', 'y5'}
    print(f"Let's test the committee W_min = {sorted(list(W_min))}")

    # Check constraints
    check1_min = len(W_min.intersection(U_X))
    check2_min = len(W_min.intersection(U_Y))
    print(f"Check 1: |W_min ∩ U_X| = {check1_min} (must be >= {l_X}) -> {'OK' if check1_min >= l_X else 'FAIL'}")
    print(f"Check 2: |W_min ∩ U_Y| = {check2_min} (must be >= {l_Y}) -> {'OK' if check2_min >= l_Y else 'FAIL'}")
    
    # Calculate representation for Voter 1
    min_approved = len(W_min.intersection(C_V1))
    print(f"Number of candidates approved by Voter 1: |W_min ∩ C_V1| = {len(W_min.intersection(C_V1))}")
    print(f"Final Minimum: {min_approved}\n")
    
    # --- 6. Final Result ---
    print("--- Conclusion ---")
    print(f"The minimum number of candidates approved by Voter 1 is {min_approved}.")
    print(f"The maximum number of candidates approved by Voter 1 is {max_approved}.")


if __name__ == '__main__':
    solve_ejr_committee_problem()