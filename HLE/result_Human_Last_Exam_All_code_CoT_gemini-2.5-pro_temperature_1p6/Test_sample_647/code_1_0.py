import math

def solve_election_ratio():
    """
    Solves the committee election problem by calculating the ratio of satisfactions.
    """
    # Voter and committee parameters
    n = 10  # Total number of voters
    k = 20  # Committee size
    group_N_size = 8 # Number of voters in group N

    # --- Step 1: Calculate s(N, W1), the lowest satisfaction for a core committee ---

    print("Step 1: Calculating the minimum satisfaction for a Core committee, s(N, W1).")
    # A property of a core committee is Proportionality for Solid Coalitions (PSC).
    # Group N = {1,...,8} is a solid coalition. They all approve candidates C_common = {1,...,8}.
    # Let's check the PSC condition for N: |N| > (l-1)*n/k
    # Here, |N|=8, n=10, k=20, so n/k = 0.5. Let's test for l = 8 (for the 8 common candidates).
    # 8 > (8-1) * 0.5  =>  8 > 3.5. This is true.
    # PSC guarantees that a core committee must contain at least l=8 candidates from C_common.
    # Since |C_common| = 8, this means any core committee W must contain all of {1,...,8}.
    print("Based on the Proportionality for Solid Coalitions (PSC) property of the core,")
    print("any core committee W1 must contain all 8 common candidates {1, ..., 8}.")

    # The base satisfaction from these 8 candidates for the 8 voters in N is:
    # Each of the 8 voters approves these 8 candidates.
    s_n_w1_base = group_N_size * 8
    print(f"The baseline satisfaction from these 8 members is {group_N_size} * 8 = {s_n_w1_base}.")

    # Now we need to choose the remaining k - 8 = 12 members of the committee.
    # To minimize the satisfaction for N, we should pick candidates that members of N do not approve.
    # Candidates {25, ..., 32} are approved by voters 9 and 10, but not by any voter in N. There are 8 such candidates.
    # We add these 8 candidates to W1. This does not increase s(N, W1).
    remaining_slots = k - 8
    non_n_candidates = 8
    remaining_slots -= non_n_candidates
    print(f"To minimize satisfaction, we select {non_n_candidates} candidates ({25,...,32}) not approved by group N.")
    
    # We still need to fill 12 - 8 = 4 slots.
    # These must be chosen from the remaining candidates, {9, ..., 24}.
    # Each candidate in {9, ..., 24} is approved by exactly one voter in N.
    # Therefore, choosing any 4 candidates from this set will add 4 to the total satisfaction score.
    added_satisfaction_w1 = remaining_slots
    print(f"The final {remaining_slots} slots must be filled by candidates from {{9,...,24}}.")
    print(f"Each of these adds 1 to the total satisfaction, so we add {added_satisfaction_w1}.")

    s_n_w1 = s_n_w1_base + added_satisfaction_w1
    print(f"So, s(N, W1) = {s_n_w1_base} + {added_satisfaction_w1} = {s_n_w1}.\n")

    # --- Step 2: Calculate s(N, W2), the lowest satisfaction for an EJR committee ---

    print("Step 2: Calculating the minimum satisfaction for an EJR committee, s(N, W2).")
    # For Extended Justified Representation (EJR), consider the group N and l=8 candidates {1,...,8}.
    # Since |N| = 8 >= n/k = 0.5, EJR applies.
    # EJR states that it cannot be that W is disjoint from {1..8} AND |A(i) intersect W| < 8 for all i in N.
    # If W is disjoint from {1..8}, then for any i in N, |A(i) intersect W| can be at most 2.
    # Since 2 < 8, this means that to satisfy EJR, W must NOT be disjoint from {1,...,8}.
    # So, any EJR committee must contain at least one candidate from {1,...,8}.
    print("Based on the Extended Justified Representation (EJR) property,")
    print("any EJR committee W2 must contain at least one candidate from {1, ..., 8}.")

    # To minimize satisfaction for N, we assume W2 contains exactly one candidate from {1,...,8}.
    # This one candidate is approved by all 8 voters in N.
    # The base satisfaction is:
    s_n_w2_base = group_N_size * 1
    print(f"To minimize satisfaction, we include just 1 of these candidates.")
    print(f"The baseline satisfaction is {group_N_size} * 1 = {s_n_w2_base}.")
    
    # We need to choose the remaining k - 1 = 19 members.
    # Again, we first pick the 8 candidates from {25, ..., 32} not approved by N.
    remaining_slots = k - 1
    remaining_slots -= non_n_candidates
    print(f"We select {non_n_candidates} candidates ({25,...,32}) not approved by group N.")
    
    # We still need to fill 19 - 8 = 11 slots.
    # These must be chosen from {9, ..., 24} to minimize added satisfaction.
    # Each adds 1 to the total satisfaction.
    added_satisfaction_w2 = remaining_slots
    print(f"The final {remaining_slots} slots are filled by candidates from {{9,...,24}}, adding {added_satisfaction_w2} to the satisfaction.")
    
    s_n_w2 = s_n_w2_base + added_satisfaction_w2
    print(f"So, s(N, W2) = {s_n_w2_base} + {added_satisfaction_w2} = {s_n_w2}.\n")

    # --- Step 3: Compute the ratio ---
    print("Step 3: Calculating the final ratio.")
    if s_n_w2 == 0:
        ratio = float('inf')
        print("s(N, W2) is zero, ratio is infinite.")
    else:
        ratio = s_n_w1 / s_n_w2
        print(f"The ratio s(N,W1) / s(N,W2) is {s_n_w1} / {s_n_w2} = {ratio:.4f}")
    
    return ratio

# Run the solver and print the final answer in the required format.
final_ratio = solve_election_ratio()
print(f"\n<<<68/19>>>")
