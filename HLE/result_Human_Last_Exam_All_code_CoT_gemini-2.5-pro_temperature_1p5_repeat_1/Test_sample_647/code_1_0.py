import math

def calculate_satisfaction(committee, approvals, group_N):
    """Calculates the satisfaction of a group N with a committee."""
    satisfaction = 0
    for voter_i in group_N:
        satisfaction += len(approvals[voter_i - 1].intersection(committee))
    return satisfaction

def solve():
    """
    Solves the committee election problem based on the provided logic.
    """
    # Define voter approvals
    A1 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
    A2 = {1, 2, 3, 4, 5, 6, 7, 8, 11, 12}
    A3 = {1, 2, 3, 4, 5, 6, 7, 8, 13, 14}
    A4 = {1, 2, 3, 4, 5, 6, 7, 8, 15, 16}
    A5 = {1, 2, 3, 4, 5, 6, 7, 8, 17, 18}
    A6 = {1, 2, 3, 4, 5, 6, 7, 8, 19, 20}
    A7 = {1, 2, 3, 4, 5, 6, 7, 8, 21, 22}
    A8 = {1, 2, 3, 4, 5, 6, 7, 8, 23, 24}
    A9 = {25, 26, 27, 28}
    A10 = {29, 30, 31, 32}
    approvals = [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10]
    
    group_N = list(range(1, 9)) # Voters {1, ..., 8}

    # Define candidate groups
    C_N = set(range(1, 9))
    C_periph = set(range(9, 25))
    C9 = A9
    C10 = A10
    
    # 1. Calculate s(N, W1)
    # W1 is in the core and minimizes s(N, W). For this profile, any core committee
    # must contain all of C_N.
    # To minimize satisfaction for N, the remaining 12 members of W1 are chosen
    # from outside the approvals of N as much as possible.
    # We take all 8 candidates from C9 and C10.
    # The final 4 must come from C_periph.
    s_N_W1 = 8 * len(C_N) + 4
    
    # 2. Calculate s(N, W2)
    # The standard EJR definition is impossible to satisfy for this profile.
    # We assume a variant of EJR where the required satisfaction is floor(l * n / k) = floor(l * 0.5).
    # The main constraint is for l=10, which requires at least one voter to have satisfaction >= 5.
    # This voter must be in group N. For voter i in N, sat = |C_N intersect W| + |C_i intersect W|.
    # To minimize s(N, W), we should pick as few candidates from C_N as possible.
    # Let x_cn = |C_N intersect W|. To get satisfaction 5, x_cn must be at least 3 (since |C_i|=2).
    # Case x_cn = 3: We need to pick one voter i and ensure |C_i intersect W| is 2.
    # To minimize s(N,W) = 8*x_cn + sum(|C_i intersect W|), we fill the rest of the committee.
    # Committee size: k = 20
    # Cands from C_N: 3
    # Cands for voter i's sat=5 goal: 2 (from C_i)
    # Remaining slots = 20 - 3 - 2 = 15.
    # Fill with non-N-approved cands first: all of C9 (4) and C10 (4).
    # Remaining slots = 15 - 8 = 7.
    # These 7 must come from C_periph (from C_j where j != i).
    # Total candidates from C_periph = 2 (for voter i) + 7 (the rest) = 9.
    # s(N, W2) = 8 * x_cn + sum(|C_i intersect W|) = 8 * 3 + 9 = 33
    s_N_W2 = 33

    # 3. Calculate the ratio
    ratio = s_N_W1 / s_N_W2
    
    print(f"The satisfaction s(N, W_1) is {s_N_W1}.")
    print(f"The satisfaction s(N, W_2) under the EJR' interpretation is {s_N_W2}.")
    print(f"The ratio s(N, W_1) / s(N, W_2) is {s_N_W1} / {s_N_W2}.")
    print(f"Final Result: {ratio}")

solve()