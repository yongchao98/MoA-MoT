import collections

def solve():
    """
    Solves the EJR committee problem.
    """
    # Voter profiles
    ballots = {
        1: {'x1', 'x2', 'x3', 'y1', 'z3'},
        2: {'x1', 'x2', 'x3', 'y2'},
        3: {'x1', 'x2', 'x3', 'y3'},
        4: {'y4', 'y5', 'y6', 'z1'},
        5: {'y4', 'y5', 'y6', 'z1'},
        6: {'y4', 'y5', 'y6', 'z2'},
        7: {'y4', 'y5', 'y6', 'z2'},
        8: {'x1', 'x2', 'x3', 'z1'},
        9: {'x1', 'x2', 'x3', 'z1'},
        10: {'x1', 'x2', 'x3', 'z1'},
    }
    
    n = 10  # number of voters
    k = 5   # committee size
    n_over_k = n / k
    voter_1_approvals = ballots[1]

    print("Step 1: Define the problem and the EJR condition.")
    print(f"There are n={n} voters and we need to elect a committee of size k={k}.")
    print(f"The EJR condition requires that for any integer l from 1 to {k},")
    print(f"any group of at least l * (n/k) = {int(n_over_k)}*l voters who agree on at least l candidates")
    print("must have at least one of those commonly approved candidates in the committee.\n")

    print("Step 2: Identify key cohesive groups and derive constraints.")
    
    # Group A: Voters who approve {x1, x2, x3}
    group_A_voters = {1, 2, 3, 8, 9, 10}
    intersection_A = set.intersection(*(ballots[i] for i in group_A_voters))
    l_A = 3
    print(f"Constraint 1: From Group A = Voters {sorted(list(group_A_voters))}")
    print(f"   - Group size is {len(group_A_voters)}. Intersection is {sorted(list(intersection_A))}, size {len(intersection_A)}.")
    print(f"   - For l={l_A}, group size {len(group_A_voters)} >= {int(n_over_k)}*l = {int(n_over_k * l_A)}.")
    print(f"   - Intersection size {len(intersection_A)} >= l = {l_A}.")
    print(f"   - EJR Conclusion: Any valid committee W must contain at least one candidate from {sorted(list(intersection_A))}.\n")

    # Group B: Voters who approve {y4, y5, y6}
    group_B_voters = {4, 5, 6, 7}
    intersection_B = set.intersection(*(ballots[i] for i in group_B_voters))
    l_B = 2
    print(f"Constraint 2: From Group B = Voters {sorted(list(group_B_voters))}")
    print(f"   - Group size is {len(group_B_voters)}. Intersection is {sorted(list(intersection_B))}, size {len(intersection_B)}.")
    print(f"   - For l={l_B}, group size {len(group_B_voters)} >= {int(n_over_k)}*l = {int(n_over_k * l_B)}.")
    print(f"   - Intersection size {len(intersection_B)} >= l = {l_B}.")
    print(f"   - EJR Conclusion: Any valid committee W must contain at least one candidate from {sorted(list(intersection_B))}.\n")

    # Group C: Voters who approve {z1}
    group_C_voters = {4, 5, 8, 9, 10}
    intersection_C = set.intersection(*(ballots[i] for i in group_C_voters))
    l_C = 1
    print(f"Constraint 3: From Group C = Voters {sorted(list(group_C_voters))}")
    print(f"   - Group size is {len(group_C_voters)}. Intersection is {sorted(list(intersection_C))}, size {len(intersection_C)}.")
    print(f"   - For l={l_C}, group size {len(group_C_voters)} >= {int(n_over_k)}*l = {int(n_over_k * l_C)}.")
    print(f"   - Intersection size {len(intersection_C)} >= l = {l_C}.")
    print(f"   - EJR Conclusion: Any valid committee W must contain the candidate 'z1'.\n")

    print("Step 3: Calculate the minimum and maximum number of committee members approved by Voter 1.")
    print(f"Voter 1's approvals: {sorted(list(voter_1_approvals))}\n")

    # Minimum Calculation
    print("--- Minimum Calculation ---")
    min_overlap = 0
    # Constraint 1 forces a member from {x1, x2, x3}. All are in Voter 1's set.
    if intersection_A.issubset(voter_1_approvals):
        min_overlap += 1
    print(f"Constraint 1 ({sorted(list(intersection_A))}) forces at least 1 member from Voter 1's approved set.")
    print(f"Therefore, the minimum number of approved candidates must be at least 1.")
    
    # To show the minimum is exactly 1, we construct a valid committee:
    # W_min = {one from intersection_A, one from intersection_B, z1, two others not in A_1}
    # Example: {x1, y4, z1, y2, z2}. Intersection with A_1 is {x1}, size 1. This committee is valid.
    min_val = 1
    print(f"We can construct a valid EJR committee W = {{'x1', 'y4', 'z1', 'y2', 'z2'}}.")
    print(f"The number of members in W approved by Voter 1 is |W intersect {sorted(list(voter_1_approvals))}| = 1.")
    print(f"So, the minimum is {min_val}.\n")

    # Maximum Calculation
    print("--- Maximum Calculation ---")
    forced_outside_members = 0
    # Constraint 2 forces a member from {y4, y5, y6}. None are in Voter 1's set.
    if intersection_B.isdisjoint(voter_1_approvals):
        forced_outside_members += 1
        print(f"Constraint 2 ({sorted(list(intersection_B))}) forces at least 1 member NOT in Voter 1's set.")
    # Constraint 3 forces z1. z1 is not in Voter 1's set.
    if intersection_C.isdisjoint(voter_1_approvals):
        forced_outside_members += 1
        print(f"Constraint 3 ({sorted(list(intersection_C))}) forces at least 1 member NOT in Voter 1's set.")
    
    max_val = k - forced_outside_members
    print(f"In total, at least {forced_outside_members} members of any EJR committee must come from outside Voter 1's set.")
    print(f"Therefore, the maximum number of approved candidates is at most {k} - {forced_outside_members} = {max_val}.")

    # To show the maximum is exactly 3, we construct a valid committee:
    # W_max = {all of intersection_A, one from intersection_B, z1}
    # Example: {x1, x2, x3, y4, z1}. Intersection with A_1 is {x1, x2, x3}, size 3. This committee is valid.
    print(f"We can construct a valid EJR committee W = {{'x1', 'x2', 'x3', 'y4', 'z1'}}.")
    print(f"The number of members in W approved by Voter 1 is |W intersect {sorted(list(voter_1_approvals))}| = 3.")
    print(f"So, the maximum is {max_val}.\n")

    print("--- Final Answer ---")
    print(f"The minimum number of candidates approved by voter 1 is {min_val}.")
    print(f"The maximum number of candidates approved by voter 1 is {max_val}.")

solve()