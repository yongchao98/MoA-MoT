import itertools

def solve_ejr_problem():
    """
    Solves the EJR problem by analyzing cohesive groups and their constraints.
    """
    # Step 0: Define the problem setup
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
    n = len(voters)
    k = 5
    voter1_approval_set = voters[1]

    print("Step 1: Understand EJR parameters")
    print(f"Number of voters (n) = {n}")
    print(f"Committee size (k) = {k}")
    print(f"The ratio n/k is {n}/{k} = {n/k}")
    print("A group V' is l-cohesive if |V'| >= l * (n/k) and its members unanimously approve >= l candidates.")
    print("An EJR committee W must satisfy |W ∩ (union of approvals of V')| >= l for every l-cohesive group V'.\n")

    print("Step 2: Identify key cohesive groups and their EJR constraints")
    # Group G1
    g1_indices = {1, 2, 3, 8, 9, 10}
    g1_ballots = [voters[i] for i in g1_indices]
    g1_intersection = set.intersection(*g1_ballots)
    g1_union = set.union(*g1_ballots)
    l_g1 = 3
    print(f"Group G1 = V{sorted(list(g1_indices))} is a {l_g1}-cohesive group:")
    print(f"  - Size condition: |G1|={len(g1_indices)} is >= {l_g1} * {n/k} = {l_g1 * (n/k)} -> True")
    print(f"  - Intersection condition: |∩G1|={len(g1_intersection)} is >= {l_g1} -> True")
    print(f"EJR Constraint 1: Any EJR committee W must select at least {l_g1} candidates from G1's union of approvals: {g1_union}\n")

    # Group G2
    g2_indices = {4, 5, 6, 7}
    g2_ballots = [voters[i] for i in g2_indices]
    g2_intersection = set.intersection(*g2_ballots)
    g2_union = set.union(*g2_ballots)
    l_g2 = 2
    print(f"Group G2 = V{sorted(list(g2_indices))} is a {l_g2}-cohesive group:")
    print(f"  - Size condition: |G2|={len(g2_indices)} is >= {l_g2} * {n/k} = {l_g2 * (n/k)} -> True")
    print(f"  - Intersection condition: |∩G2|={len(g2_intersection)} is >= {l_g2} -> True")
    print(f"EJR Constraint 2: Any EJR committee W must select at least {l_g2} candidates from G2's union of approvals: {g2_union}\n")

    print("Step 3: Calculate the Maximum number of approved candidates for Voter 1")
    voter1_g2_union_intersection = voter1_approval_set.intersection(g2_union)
    print(f"Voter 1's approval set A1 = {voter1_approval_set}")
    print(f"Intersection of A1 and G2's approval union: A1 ∩ U(G2) = {voter1_g2_union_intersection}, size = {len(voter1_g2_union_intersection)}")
    print(f"To satisfy Constraint 2, W must have at least {l_g2} members from {g2_union}.")
    print("Since these candidates are not approved by Voter 1, this leaves at most k - l_g2 seats for other candidates.")
    max_val = k - l_g2
    print(f"Max candidates for Voter 1 <= {k} - {l_g2} = {max_val}")
    print("Let's test if a committee can achieve this maximum of 3.")
    W_max = {'x1', 'x2', 'x3', 'y4', 'y5'}
    print(f"Example committee W_max = {W_max}")
    print(f"Check C1: |W_max ∩ U(G1)| = {len(W_max.intersection(g1_union))} >= {l_g1} -> True")
    print(f"Check C2: |W_max ∩ U(G2)| = {len(W_max.intersection(g2_union))} >= {l_g2} -> True")
    print(f"Intersection with Voter 1's approvals: |W_max ∩ A1| = {len(W_max.intersection(voter1_approval_set))}")
    print(f"The maximum is indeed {max_val}.\n")

    print("Step 4: Calculate the Minimum number of approved candidates for Voter 1")
    print("Let's try to construct an EJR committee with 0 members approved by Voter 1.")
    W_min = {'z1', 'y2', 'y3', 'y4', 'y5'}
    print(f"Proposed committee W_min = {W_min}")
    print(f"Check C1: |W_min ∩ U(G1)| = {len(W_min.intersection(g1_union))} >= {l_g1} -> True")
    print(f"Check C2: |W_min ∩ U(G2)| = {len(W_min.intersection(g2_union))} >= {l_g2} -> True")
    min_val_intersect_size = len(W_min.intersection(voter1_approval_set))
    print(f"Intersection with Voter 1's approvals: |W_min ∩ A1| = {min_val_intersect_size}")
    min_val = min_val_intersect_size
    print(f"The minimum is {min_val}.\n")

    print(f"Conclusion: The minimum number of candidates approved by Voter 1 is {min_val} and the maximum is {max_val}.")

solve_ejr_problem()
<<<L>>>