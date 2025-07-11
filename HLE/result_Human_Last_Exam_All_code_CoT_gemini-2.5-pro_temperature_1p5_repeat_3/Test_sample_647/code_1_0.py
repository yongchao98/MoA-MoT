import math

def solve():
    """
    Solves the committee election problem by calculating the satisfaction
    for committees W1 (core) and W2 (EJR) and finding their ratio.
    """

    # --- Step 1: Define the problem setup ---
    candidates = set(range(1, 33))
    C_core = set(range(1, 9))
    C_N_spec = set(range(9, 25))
    C_others = set(range(25, 33))
    
    # N is the set of indices for voters 1 to 8
    voter_indices_N = list(range(8)) # 0-indexed for a list

    committee_size = 20
    total_voters = 10

    # --- Step 2: Calculate s(N, W1) for a committee in the Core ---
    # A committee in the core maximizes the total approval score.
    # Approval scores:
    # - Candidates in C_core ({1..8}) are approved by 8 voters. Score = 8.
    # - All other candidates are approved by 1 voter. Score = 1.
    # To maximize the score, any core committee of size 20 must contain the 8
    # candidates from C_core. The other 12 members are chosen from the rest.
    
    W1_base = C_core
    s_N_W1_base = 8 * len(W1_base) # Each of the 8 voters in N approves all 8 members

    # To minimize s(N, W), we must choose the remaining 12 members of the
    # committee (W_rem) to minimize the additional satisfaction.
    # The additional satisfaction is sum(|A(i) intersect W_rem|) for i in N.
    # This is equivalent to minimizing |C_N_spec intersect W_rem|.
    
    # Candidates for W_rem are those not in C_core.
    rem_candidates = C_N_spec.union(C_others)
    
    # To minimize intersection with C_N_spec, we should prioritize candidates
    # from C_others.
    num_from_C_others = len(C_others) # which is 8
    
    # We must pick 12 members for the remainder of the committee.
    # We pick all 8 from C_others.
    # The remaining 12 - 8 = 4 must be chosen from C_N_spec.
    num_from_C_N_spec = 4
    
    # These 4 candidates from C_N_spec will each add 1 to the satisfaction sum.
    s_N_W1_additional = num_from_C_N_spec
    
    s_N_W1 = s_N_W1_base + s_N_W1_additional
    
    # --- Step 3: Calculate s(N, W2) for a committee satisfying EJR ---
    # We analyze the constraints imposed by EJR.
    # n/k = 10/20 = 0.5
    #
    # Constraint 1 (Group N={1..8}): |N|=8. They are 8-cohesive (|C_core|=8).
    # Since |N| = 8 >= l * (n/k) for l=8 (8 >= 8*0.5=4), EJR requires
    # max_{i in N} |A(i) intersect W| >= 8.
    #
    # Constraint 2 (Voter 9): |{9}|=1. Voter 9 is 2-cohesive (|A(9)|=4).
    # Since |{9}|=1 >= l*(n/k) for l=2 (1 >= 2*0.5=1), EJR requires
    # |A(9) intersect W| >= 2. So W must have at least 2 members from C_others_9={25..28}.
    #
    # Constraint 3 (Voter 10): Similarly, EJR requires |A(10) intersect W| >= 2.
    # W must have at least 2 members from C_others_10={29..32}.
    #
    # This implies a committee W must contain at least 4 members from C_others.
    
    # We want to minimize s(N, W) = 8 * |W intersect C_core| + |W intersect C_N_spec|.
    # Let x = |W intersect C_core| and y = |W intersect C_N_spec|. Minimize 8x + y.
    
    min_s_N_W2 = float('inf')
    
    # Iterate through possible numbers of candidates from C_core (x)
    for x in range(len(C_core) + 1):
        # Check the EJR constraint for group N.
        # max |A(i) intersect W| = x + max |C_i' intersect W_n| >= 8.
        # Since |C_i'|=2, max |C_i' intersect W_n| can be at most 2.
        # So x + 2 >= 8, which means x must be at least 6.
        if x < 6:
            continue
            
        # Determine the minimum y required for the EJR constraint for N
        if x == 7: # Need max |C_i' intersect W_n| >= 1, so we need y >= 1.
            min_y_for_ejr_N = 1
        elif x == 6: # Need max |C_i' intersect W_n| >= 2, so we need y >= 2.
            min_y_for_ejr_N = 2
        else: # x == 8
            min_y_for_ejr_N = 0
            
        # Iterate through possible numbers of candidates from C_N_spec (y)
        for y in range(len(C_N_spec) + 1):
            if y < min_y_for_ejr_N:
                continue

            # Number of candidates from C_others (z) is determined by k.
            z = committee_size - x - y
            
            # Check constraints on z
            # 1. z must be non-negative and at most len(C_others)
            if not (0 <= z <= len(C_others)):
                continue
            
            # 2. z must satisfy EJR for voters 9 and 10, so z >= 4
            if z < 4:
                continue
            
            # This combination (x, y, z) is valid. Calculate satisfaction.
            current_s = 8 * x + y
            if current_s < min_s_N_W2:
                min_s_N_W2 = current_s
                
    s_N_W2 = min_s_N_W2
    
    # --- Step 4: Final calculation and output ---
    print(f"The lowest satisfaction for N with a core committee is s(N, W1) = {s_N_W1}")
    print(f"The lowest satisfaction for N with an EJR committee is s(N, W2) = {s_N_W2}")
    
    # Use math.gcd to simplify the fraction
    common_divisor = math.gcd(s_N_W1, s_N_W2)
    s1_simple = s_N_W1 // common_divisor
    s2_simple = s_N_W2 // common_divisor

    print(f"The ratio is s(N, W1) / s(N, W2) = {s_N_W1} / {s_N_W2} = {s1_simple} / {s2_simple}")

solve()
<<<1.2592592592592593>>>