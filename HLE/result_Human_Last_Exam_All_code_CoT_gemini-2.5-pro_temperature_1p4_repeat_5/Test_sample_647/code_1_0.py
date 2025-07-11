import math

def solve():
    """
    Solves the committee election problem by calculating the satisfaction scores
    for the two committees W1 and W2 and then finding their ratio.
    """

    # --- Step 1: Define voters' approval sets and group N ---
    A = {
        1: {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
        2: {1, 2, 3, 4, 5, 6, 7, 8, 11, 12},
        3: {1, 2, 3, 4, 5, 6, 7, 8, 13, 14},
        4: {1, 2, 3, 4, 5, 6, 7, 8, 15, 16},
        5: {1, 2, 3, 4, 5, 6, 7, 8, 17, 18},
        6: {1, 2, 3, 4, 5, 6, 7, 8, 19, 20},
        7: {1, 2, 3, 4, 5, 6, 7, 8, 21, 22},
        8: {1, 2, 3, 4, 5, 6, 7, 8, 23, 24},
        9: {25, 26, 27, 28},
        10: {29, 30, 31, 32},
    }
    N = list(range(1, 9))
    committee_size = 20

    # --- Step 2: Analyze the Core and find s(N, W1) ---
    # A committee W is in the core if it's not "blocked". It is blocked by a pair (N', W')
    # if W' is a committee of size 20, W' is a subset of candidates approved by N',
    # and N' gets higher satisfaction from W' than from W.

    # Let's consider the coalition N' = N = {1, ..., 8}.
    # The set of candidates they can form a committee from is union_{i in N} A(i) = {1, ..., 24}.
    # To maximize their satisfaction, they would form a committee W_best_for_N by picking
    # the 8 candidates common to all of them, {1,...,8}, and 12 other candidates from {9,...,24}.
    # The satisfaction of N with such a committee would be:
    # Each of the 8 members of N approves the 8 common candidates, contributing 8*8=64.
    # Each of the 12 other candidates is approved by exactly one member of N, contributing 12*1=12.
    # Max satisfaction for N is 64 + 12 = 76.
    max_satisfaction_for_N = 8 * 8 + 12 * 1
    
    # Any committee W for which s(N, W) < 76 is blocked by N. For example, if W gave N
    # a satisfaction of 75, N could propose W_best_for_N, from which they get satisfaction 76,
    # and block W.
    # Therefore, any committee W1 in the core must satisfy s(N, W1) >= 76.
    # It can be shown that committees achieving this value (e.g., W={1,...,20}) are not blocked
    # by other coalitions and are thus in the core. This makes the satisfaction for N of any core committee
    # equal to 76.
    s_N_W1 = max_satisfaction_for_N
    
    print(f"Step 1: The lowest satisfaction for group N from a core committee, s(N,W1).")
    print(f"Any core committee must provide group N a satisfaction of at least {s_N_W1}.")
    print(f"Thus, s(N,W1) = {s_N_W1}")
    print("-" * 20)

    # --- Step 3: Analyze EJR and find s(N, W2) ---
    # EJR requires that for any group of l voters N', if |intersection(A(i) for i in N')| >= l * k/n,
    # the committee must select at least one candidate from that intersection.
    # k=20, n=10 => k/n = 2.
    
    # EJR constraints:
    # 1. For N'={1,2,3,4} (l=4), intersection is {1..8}. |{1..8}|=8. Since 8 >= 4*2, W must have a candidate from {1..8}.
    # 2. For N'={9} (l=1), intersection is A(9). |A(9)|=4. Since 4 >= 1*2, W must have a candidate from {25..28}.
    # 3. For N'={10} (l=1), intersection is A(10). |A(10)|=4. Since 4 >= 1*2, W must have a candidate from {29..32}.
    
    # To find W2 that minimizes s(N,W), we build a committee satisfying these rules with candidates
    # that are least approved by N.
    # Contribution to s(N,W) by candidate c: 8 if c in {1..8}, 1 if c in {9..24}, 0 if c in {25..32}.
    
    s_N_W2 = 0
    
    # From constraint 1, we must pick one candidate from {1..8}. To minimize satisfaction, we pick one.
    # This candidate adds 8 to the satisfaction score for group N.
    satisfaction_from_common = 8
    s_N_W2 += satisfaction_from_common
    
    # From constraints 2 and 3, we must pick candidates from {25..32}. These add 0 to satisfaction.
    # We can pick all 8 of them ({25..32}) to fill committee slots with "zero-cost" candidates.
    satisfaction_from_v9_v10_cands = 0
    s_N_W2 += satisfaction_from_v9_v10_cands

    # We have now chosen 1 (from {1..8}) + 8 (from {25..32}) = 9 candidates.
    # We need to choose 20 - 9 = 11 more.
    # To continue minimizing satisfaction, we choose from the next cheapest pool: {9..24}, where each adds 1.
    satisfaction_from_pairs = 11 * 1
    s_N_W2 += satisfaction_from_pairs
    
    print(f"Step 2: The lowest satisfaction for group N from an EJR committee, s(N,W2).")
    print(f"An EJR committee W2 must contain at least one candidate from {{1..8}}, which contributes {satisfaction_from_common} to the satisfaction of N.")
    print(f"To minimize s(N,W2), the other 19 members are chosen from candidates least approved by N.")
    print(f"The committee will contain all 8 candidates from {{25..32}} (0 satisfaction) and 11 from {{9..24}} ({satisfaction_from_pairs} satisfaction).")
    print(f"Total satisfaction s(N,W2) = {satisfaction_from_common} + {satisfaction_from_v9_v10_cands} + {satisfaction_from_pairs} = {s_N_W2}")
    print("-" * 20)

    # --- Step 4: Compute the Final Ratio ---
    ratio = s_N_W1 / s_N_W2
    print(f"Step 3: Calculating the ratio s(N,W1)/s(N,W2).")
    print(f"The ratio is {s_N_W1} / {s_N_W2} = {ratio}")
    
    return ratio

if __name__ == '__main__':
    final_ratio = solve()
    # The final answer is wrapped in <<<>>>
    print(f"\n<<<{final_ratio}>>>")
