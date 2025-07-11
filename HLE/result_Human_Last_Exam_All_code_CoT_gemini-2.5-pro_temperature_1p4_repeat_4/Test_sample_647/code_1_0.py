import math

def solve_election_problem():
    """
    Solves the committee election problem by calculating the satisfactions s(N, W1) and s(N, W2)
    and then their ratio.
    """

    # Define candidate sets
    C_shared = set(range(1, 9))  # {1, ..., 8}, approved by voters 1-8
    C_pairs = {i: set(range(9 + (i-1)*2, 9 + (i-1)*2 + 2)) for i in range(1, 9)} # {9,10}, {11,12}, ...
    C_lone_1 = set(range(25, 29)) # {25,26,27,28}, approved by voter 9
    C_lone_2 = set(range(29, 33)) # {29,30,31,32}, approved by voter 10
    
    # Define voter approval sets
    A = {i: C_shared.union(C_pairs[i]) for i in range(1, 9)}
    A[9] = C_lone_1
    A[10] = C_lone_2
    
    group_N = range(1, 9)

    # Part 1: Find s(N, W1)
    # W1 is a committee of size 20 in the core that minimizes satisfaction for group N.
    # The core contains AV winners. The 8 candidates in C_shared have an approval score of 8.
    # All other 24 candidates have a score of 1.
    # Thus, any core committee must contain C_shared.
    # W1 = C_shared U C_other, where |C_other| = 12.
    # To minimize s(N, W1), we must choose C_other to have the minimum overlap with
    # the candidates approved by N, which are C_shared U (union of C_pairs).
    # This means we should pick candidates from C_lone_1 and C_lone_2 first.
    # C_other should contain all of C_lone_1 (4) and C_lone_2 (4).
    # We need to pick 4 more candidates. They must come from C_pairs.
    # Let's pick them from C_pairs[1] and C_pairs[2].
    
    W1 = C_shared.union(C_lone_1).union(C_lone_2).union(C_pairs[1]).union(C_pairs[2])

    s_N_W1_list = []
    for i in group_N:
        s_N_W1_list.append(len(A[i].intersection(W1)))
    
    s_N_W1 = sum(s_N_W1_list)

    print("--- Calculation for s(N, W1) ---")
    print(f"The committee W1 which is in the core and minimizes satisfaction for N is composed of:")
    print(f"- The 8 shared candidates: {sorted(list(C_shared))}")
    print(f"- The 8 candidates from voters 9 and 10: {sorted(list(C_lone_1.union(C_lone_2)))}")
    print(f"- 4 candidates from the pairs to fill the committee: {sorted(list(C_pairs[1].union(C_pairs[2])))}")
    print("Satisfaction for each voter in N = {1..8} with W1:")
    print(f"s(N, W1) = {' + '.join(map(str, s_N_W1_list))} = {s_N_W1}")
    print("-" * 35)
    
    # Part 2: Find s(N, W2)
    # W2 is a committee of size 20 satisfying EJR that minimizes satisfaction for group N.
    # EJR constraints:
    # 1. For group {9} (l=1, cohesive), |A(9) intersect W2| >= 1.
    # 2. For group {10} (l=1, cohesive), |A(10) intersect W2| >= 1.
    # 3. For group N={1..8} (l=8, cohesive), exists i in N with |A(i) intersect W2| >= 8.
    # To minimize s(N,W2), we should first pick all candidates from C_lone_1 and C_lone_2 (8 candidates).
    # This satisfies constraints 1 and 2.
    # We need to pick 12 more candidates from C_shared U C_pairs. Let this set be W_prime.
    # We must satisfy constraint 3: exists i, |A(i) intersect W_prime| >= 8.
    # |A(i) intersect W_prime| = |C_shared intersect W_prime| + |C_pair_i intersect W_prime|.
    # To minimize overall satisfaction for N, we must select as few candidates from C_shared as possible.
    # The minimum number of shared candidates is 6. This is achieved by taking 2 candidates
    # from a C_pair_i (e.g., C_pairs[1]) and 6 from C_shared.
    # This selection satisfies the EJR constraint for N.
    # The remaining 12 - 8 = 4 candidates for W_prime should be from other C_pairs to keep satisfaction low.
    
    W2_c_lones = C_lone_1.union(C_lone_2)
    W2_c_shared_part = set(range(1, 7)) # {1,2,3,4,5,6}
    W2_c_pairs_part = C_pairs[1].union(C_pairs[2]).union(C_pairs[3]) # {9,10,11,12,13,14}
    
    W2 = W2_c_lones.union(W2_c_shared_part).union(W2_c_pairs_part)

    s_N_W2_list = []
    for i in group_N:
        s_N_W2_list.append(len(A[i].intersection(W2)))
        
    s_N_W2 = sum(s_N_W2_list)

    print("--- Calculation for s(N, W2) ---")
    print(f"The committee W2 which satisfies EJR and minimizes satisfaction for N is composed of:")
    print(f"- The 8 candidates from voters 9 and 10: {sorted(list(W2_c_lones))}")
    print(f"- 6 candidates from C_shared to satisfy EJR for N minimally: {sorted(list(W2_c_shared_part))}")
    print(f"- 6 candidates from C_pairs to fill the committee: {sorted(list(W2_c_pairs_part))}")
    print("Satisfaction for each voter in N = {1..8} with W2:")
    print(f"s(N, W2) = {' + '.join(map(str, s_N_W2_list))} = {s_N_W2}")
    print("-" * 35)

    # Part 3: Calculate the ratio
    ratio = s_N_W1 / s_N_W2
    
    print("--- Final Ratio ---")
    print(f"The ratio s(N,W1)/s(N,W2) is {s_N_W1} / {s_N_W2}")
    print(f"The calculated fraction is 34/27")
    print(f"The final numerical result is: {ratio}")
    print("-" * 35)

    return ratio

final_ratio = solve_election_problem()
# The final answer is wrapped according to the required format.
# print(f'<<<{final_ratio}>>>')