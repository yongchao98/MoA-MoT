import collections

def solve_committee_problem():
    """
    Solves the EJR committee problem based on the provided voter profile.
    """
    # Step 1: Define the electoral setup and voter ballots
    n = 10
    k = 5
    quota = n / k

    # Ballots are 0-indexed for voters 1-10
    ballots = [
        {'x1', 'x2', 'x3', 'y1', 'z3'},  # Voter 1
        {'x1', 'x2', 'x3', 'y2'},        # Voter 2
        {'x1', 'x2', 'x3', 'y3'},        # Voter 3
        {'y4', 'y5', 'y6', 'z1'},        # Voter 4
        {'y4', 'y5', 'y6', 'z1'},        # Voter 5
        {'y4', 'y5', 'y6', 'z2'},        # Voter 6
        {'y4', 'y5', 'y6', 'z2'},        # Voter 7
        {'x1', 'x2', 'x3', 'z1'},        # Voter 8
        {'x1', 'x2', 'x3', 'z1'},        # Voter 9
        {'x1', 'x2', 'x3', 'z1'}         # Voter 10
    ]
    voter1_ballot = ballots[0]

    print("Step 1: Understanding the EJR Constraints")
    print(f"There are n={n} voters and committee size k={k}. The Hare quota is n/k = {int(quota)}.")
    print("-" * 40)

    # Step 2: Identify key cohesive groups and their constraints
    
    # Group X: Voters {1,2,3,8,9,10} -> indices {0,1,2,7,8,9}
    group_x_indices = [0, 1, 2, 7, 8, 9]
    group_x_ballots = [ballots[i] for i in group_x_indices]
    intersection_x = set.intersection(*group_x_ballots)
    union_x = set.union(*group_x_ballots)
    l_x = 3
    
    print("Step 2: Identifying Cohesive Groups")
    print("\nAnalyzing Group X (Voters 1,2,3,8,9,10):")
    print(f" - Group size = {len(group_x_indices)}. For l=3, required size is {l_x} * {int(quota)} = {l_x * int(quota)}. Condition: {len(group_x_indices)} >= {l_x * int(quota)} -> True.")
    print(f" - Common candidates = {sorted(list(intersection_x))}. Size = {len(intersection_x)}. Required: >= {l_x} -> True.")
    print(f" - Conclusion: Group X is 3-cohesive. EJR requires any committee W to satisfy: |W ∩ C_X| >= {l_x}")
    print(f"   where C_X = union of approvals = {sorted(list(union_x))}")

    # Group Y: Voters {4,5,6,7} -> indices {3,4,5,6}
    group_y_indices = [3, 4, 5, 6]
    group_y_ballots = [ballots[i] for i in group_y_indices]
    intersection_y = set.intersection(*group_y_ballots)
    union_y = set.union(*group_y_ballots)
    l_y = 2

    print("\nAnalyzing Group Y (Voters 4,5,6,7):")
    print(f" - Group size = {len(group_y_indices)}. For l=2, required size is {l_y} * {int(quota)} = {l_y * int(quota)}. Condition: {len(group_y_indices)} >= {l_y * int(quota)} -> True.")
    print(f" - Common candidates = {sorted(list(intersection_y))}. Size = {len(intersection_y)}. Required: >= {l_y} -> True.")
    print(f" - Conclusion: Group Y is 2-cohesive. EJR requires any committee W to satisfy: |W ∩ C_Y| >= {l_y}")
    print(f"   where C_Y = union of approvals = {sorted(list(union_y))}")
    print("-" * 40)

    # Step 3: Find the minimum number of approved candidates for Voter 1
    print("Step 3: Finding the Minimum for Voter 1")
    print(f"Voter 1's approvals are A1 = {sorted(list(voter1_ballot))}.")
    print("We need to build a committee W of size 5 satisfying:")
    print(f"1. |W ∩ C_X| >= {l_x}")
    print(f"2. |W ∩ C_Y| >= {l_y}")

    # To minimize overlap with A1, we pick candidates outside of A1.
    candidates_from_x_not_in_a1 = union_x - voter1_ballot
    
    # We must pick 3 members from C_X. Let's pick them from candidates_from_x_not_in_a1 if possible.
    # candidates_from_x_not_in_a1 = {'y2', 'y3', 'z1'}. This has size 3. Perfect.
    w_min_part1 = candidates_from_x_not_in_a1
    
    # Check constraint 2 for this partial committee
    # intersection with C_Y is {'z1'}. We need 2. So we need 1 more from C_Y.
    # Let's pick from C_Y, avoiding A1 (which is easy as A1 and C_Y are disjoint)
    w_min_part2 = {'y4'} # from C_Y
    
    # We have {'y2', 'y3', 'z1', 'y4'}. Need one more to make size 5. Pick another from C_Y like 'y5'.
    w_min = {'y2', 'y3', 'z1', 'y4', 'y5'}
    min_val = len(w_min.intersection(voter1_ballot))

    print("We construct W_min to avoid candidates from A1:")
    print(f"Proposed W_min = {sorted(list(w_min))}")
    print(f"Check constraints: |W_min ∩ C_X| = |{sorted(list(w_min.intersection(union_x)))}| = {len(w_min.intersection(union_x))} (>= {l_x}, OK)")
    print(f"Check constraints: |W_min ∩ C_Y| = |{sorted(list(w_min.intersection(union_y)))}| = {len(w_min.intersection(union_y))} (>= {l_y}, OK)")
    print(f"Intersection with Voter 1's ballot: |W_min ∩ A1| = |{sorted(list(w_min))} ∩ {sorted(list(voter1_ballot))}| = {min_val}")
    print(f"The minimum number is {min_val}.")
    print("-" * 40)
    
    # Step 4: Find the maximum number of approved candidates for Voter 1
    print("Step 4: Finding the Maximum for Voter 1")
    print(f"To maximize |W ∩ A1|, we must still satisfy |W ∩ C_Y| >= {l_y}.")
    a1_intersect_cy = voter1_ballot.intersection(union_y)
    print(f"Voter 1's ballot A1 has no candidates in common with C_Y (A1 ∩ C_Y = {a1_intersect_cy}).")
    print(f"So, the {l_y} members of W that must come from C_Y cannot be in A1.")
    print(f"This means at most k - l_y = {k} - {l_y} = {k - l_y} members of W can be from A1.")
    
    # Can we achieve this maximum of 3?
    # Let's pick 3 from A1, e.g. {'x1', 'x2', 'x3'}.
    # We need 2 from C_Y, e.g. {'y4', 'z1'}.
    w_max = {'x1', 'x2', 'x3', 'y4', 'z1'}
    max_val = len(w_max.intersection(voter1_ballot))

    print(f"\nWe construct W_max to include as many candidates from A1 as possible:")
    print(f"Proposed W_max = {sorted(list(w_max))}")
    print(f"Check constraints: |W_max ∩ C_X| = |{sorted(list(w_max.intersection(union_x)))}| = {len(w_max.intersection(union_x))} (>= {l_x}, OK)")
    print(f"Check constraints: |W_max ∩ C_Y| = |{sorted(list(w_max.intersection(union_y)))}| = {len(w_max.intersection(union_y))} (>= {l_y}, OK)")
    print(f"Intersection with Voter 1's ballot: |W_max ∩ A1| = |{sorted(list(w_max))} ∩ {sorted(list(voter1_ballot))}| = {max_val}")
    print(f"The maximum number is {max_val}.")
    print("-" * 40)

    print(f"Conclusion: The minimum is {min_val} and the maximum is {max_val}.")

solve_committee_problem()
<<<L>>>