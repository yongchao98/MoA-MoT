def solve_election_problem():
    """
    Solves the committee election problem by calculating s(N, W_1) and s(N, W_2).
    """
    
    # --- Step 0: Define Problem Parameters ---
    # Number of voters n=10, committee size k=20
    # Group N = {1, ..., 8}
    # C_common = {1, ..., 8}, candidates approved by all in N.
    # C_i for i in {1..8} are disjoint sets of 2 candidates each, approved only by voter i.
    # C_9 and C_10 are disjoint sets of 4 candidates each, approved by voters 9 and 10 respectively.
    
    # --- Step 1: Find s(N, W_1) for a Core committee ---
    
    # A committee is in the core if it maximizes the Proportional Approval Voting (PAV) score.
    # PAV score is maximized by greedily selecting candidates that provide the largest marginal score increase.
    # The marginal increase for adding a candidate 'c' is the sum of 1/(k_i+1) over all voters 'i' who approve 'c',
    # where k_i is the number of candidates already in the committee that voter 'i' approves.

    # 1. Initially, candidates in C_common are approved by 8 voters (group N).
    #    Their marginal score contribution is 8 * (1/1) = 8.
    #    All other candidates are approved by only one voter, giving a marginal score of 1 * (1/1) = 1.
    #    Thus, all 8 candidates from C_common are selected first.
    #    Committee so far: C_common. Size = 8.

    # 2. Now, voters in N have 8 approved candidates. The marginal gain from adding one of their
    #    private candidates (from C_1, ..., C_8) is 1/9.
    #    Voters 9 and 10 have 0 approved candidates. Adding one of their candidates gives a gain of 1.
    #    Since 1 > 1/9, candidates from C_9 and C_10 are selected next.
    #    After adding all 4 candidates from C_9 and all 4 from C_10, the committee is:
    #    W_so_far = C_common U C_9 U C_10. Size = 8 + 4 + 4 = 16.

    # 3. We need to select 4 more candidates (20 - 16 = 4).
    #    The only remaining candidates are from C_1 U ... U C_8.
    #    The marginal gain for any of these is 1/9. So, any 4 candidates from this pool can be chosen.
    
    # Let K be the set of these 4 candidates. A core committee W_1 has the form:
    # W_1 = C_common U C_9 U C_10 U K, where K is a set of 4 candidates from C_1 U ... U C_8.
    
    # The satisfaction for group N is s(N, W) = sum_{i in N} |A(i) intersect W|.
    # For any i in N, A(i) = C_common U C_i.
    # s(N, W_1) = sum_{i=1 to 8} |(C_common U C_i) intersect (C_common U C_9 U C_10 U K)|
    #           = sum_{i=1 to 8} (|C_common| + |C_i intersect K|)
    #           = 8 * |C_common| + sum_{i=1 to 8} |C_i intersect K|
    #           = 8 * 8 + |K|
    
    s_N_W1 = 8 * 8 + 4
    
    print("--- Analysis for W_1 (Core) ---")
    print(f"Any core committee W_1 must contain all 8 common candidates, and 4 additional candidates approved by voters in N.")
    print(f"Satisfaction s(N, W_1) = 8 * 8 + 4 = {s_N_W1}")
    
    # --- Step 2: Find s(N, W_2) for an EJR committee ---
    
    # EJR requires that for any group N' of size n', if they cohesively agree on l=floor(n'*k/n)
    # candidates, then at least one voter in N' must have at least l of their approved candidates elected.
    # Here n=10, k=20, so l = floor(n'*2).
    
    # EJR Constraints for this problem:
    # - For n'=1 (any voter i), l=2. |A(i)|>=2, so |A(i) intersect W| >= 2. This applies to all 10 voters.
    # - For n'=4 (any 4 voters from N), l=8. They agree on C_common, which has size 8.
    #   So, for any 4 voters in N, at least one must have |A(i) intersect W| >= 8.

    # We want to find a committee W_2 satisfying these rules that minimizes s(N,W).
    # Let w_c = |C_common intersect W_2| and n_i = |C_i intersect W_2|.
    # The goal is to minimize s(N, W_2) = 8*w_c + sum(n_i for i in {1..8}).
    # Total committee size is w_c + sum(n_i for i in {1..10}) = 20.

    # The constraint |A(i) intersect W| >= 8 for some i in any group of 4 from N implies w_c+n_i >= 8.
    # Since n_i <= 2, we must have w_c + 2 >= 8, which means w_c >= 6.
    
    # The sum of n_i for i in {1..8} must be large enough to ensure that for any 4 voters,
    # one has n_i=2 (given w_c=6). This requires at least 4 members of N to have n_i=2.
    # So, sum(n_i for i in {1..8}) must be at least 4 * 2 = 8.
    
    # The satisfaction is s(N,W) = 8*w_c + sum_{i=1..8} n_i = 8*w_c + (20 - w_c - n_9 - n_10)
    #                             = 7*w_c + 20 - n_9 - n_10.
    # To minimize this, we must minimize w_c and maximize n_9 and n_10.
    
    # The minimal w_c is 6.
    # To maximize n_9 and n_10, we must satisfy the committee size constraint while meeting the requirement on sum(n_i for i in N).
    # sum_{i=1..8} n_i >= 8  =>  20 - w_c - n_9 - n_10 >= 8  =>  w_c + n_9 + n_10 <= 12.
    # With w_c=6, we have n_9 + n_10 <= 6.
    # EJR also requires n_9 >= 2 and n_10 >= 2.
    # To maximize n_9 + n_10, we can choose n_9=3, n_10=3 (or n_9=4, n_10=2, etc.), making the sum 6.
    
    # Optimal choice for minimizing s(N, W) is:
    w_c_W2 = 6
    n9_plus_n10_W2 = 6
    sum_ni_N_W2 = 20 - w_c_W2 - n9_plus_n10_W2 # This is 20 - 6 - 6 = 8
    
    s_N_W2 = 8 * w_c_W2 + sum_ni_N_W2

    print("\n--- Analysis for W_2 (EJR) ---")
    print("To minimize satisfaction s(N,W) while satisfying EJR, we found an optimal structure:")
    print(f"Number of common candidates (w_c): {w_c_W2}")
    print(f"Number of candidates from voters {1..8} (sum(n_i)): {sum_ni_N_W2}")
    print(f"Satisfaction s(N, W_2) = 8 * {w_c_W2} + {sum_ni_N_W2} = {s_N_W2}")
    
    # --- Step 3: Calculate the Final Ratio ---
    
    ratio = s_N_W1 / s_N_W2
    
    print("\n--- Final Calculation ---")
    print(f"The ratio is s(N,W_1) / s(N,W_2)")
    print(f"The final equation is: {s_N_W1} / {s_N_W2} = {ratio}")
    
    return ratio

final_answer = solve_election_problem()
# The final result is a fraction, let's represent it with more precision.
# 68 / 56 = 17 / 14
final_answer_val = 17/14
print(f"\nFinal Answer: {final_answer_val}")

print(f"\n<<<17/14>>>")