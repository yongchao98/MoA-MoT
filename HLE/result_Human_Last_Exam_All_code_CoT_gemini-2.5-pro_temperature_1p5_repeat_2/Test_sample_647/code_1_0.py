import sys

def calculate_satisfaction(voter_profiles, group_N, committee_W):
    """Calculates the total satisfaction for a group of voters with a committee."""
    total_satisfaction = 0
    for voter_id in group_N:
        # Voter's approved candidates are at index voter_id - 1
        approved_candidates = set(voter_profiles[voter_id - 1])
        # |A(i) intersect W|
        satisfaction = len(approved_candidates.intersection(committee_W))
        total_satisfaction += satisfaction
    return total_satisfaction

def solve():
    """
    Solves the committee election problem by finding W1 and W2 and their satisfaction ratio.
    """
    # Define voter approval sets A(i)
    voter_profiles = [
        {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
        {1, 2, 3, 4, 5, 6, 7, 8, 11, 12},
        {1, 2, 3, 4, 5, 6, 7, 8, 13, 14},
        {1, 2, 3, 4, 5, 6, 7, 8, 15, 16},
        {1, 2, 3, 4, 5, 6, 7, 8, 17, 18},
        {1, 2, 3, 4, 5, 6, 7, 8, 19, 20},
        {1, 2, 3, 4, 5, 6, 7, 8, 21, 22},
        {1, 2, 3, 4, 5, 6, 7, 8, 23, 24},
        {25, 26, 27, 28},
        {29, 30, 31, 32}
    ]

    # The group of voters N
    group_N = range(1, 9)

    # Based on our analysis, we construct W1, a core committee that minimizes s(N,W).
    # W1 must contain C_N={1..8} to be in the core (unblockable by N).
    # To minimize s(N,W1), we fill the rest of the committee with non-N candidates first.
    # W1 contains: C_N (8), all of A(9) and A(10) (8), and 4 from C_P to make size 20.
    # To max out a voter in N (making it core), P1={9,10} must be in W1.
    # But this contradicts filling with non-N candidates.
    # A revised analysis for minimal core committee W1 suggests it's sufficient to max out
    # one voter's satisfaction to guarantee it is in the core.
    # W1 = C_N U P_j U C_9 U C_10 U {2 more from C_P}
    # For instance, j=1, k=2: W1 = C_N U P_1 U C_9 U C_10 U P_2. |W1| = 8+2+8+2 = 22 (too large).
    # Correct W1 construction:
    # A core committee with minimum satisfaction for N must contain all common candidates C_N (8 members).
    # To minimize satisfaction, it should contain the maximum number of candidates outside of N's approval sets.
    # These are from A(9) and A(10) (8 members).
    # The remaining 20 - 8 - 8 = 4 members must come from the private candidates of N (C_P).
    # To ensure it is in the core, at least one voter in N must have their satisfaction maxed out.
    # A(1) = {1..10}. W1 needs {1..10}. {1..8} is C_N. {9,10} is P_1.
    # W_1 requires C_N and P_1, filling 10 slots. The other 10 slots should be filled to minimize s(N, W_1),
    # so we take all 8 from C_9 and C_10, and 2 more from C_P (e.g. P_2).
    
    C_N = set(range(1, 9))
    C_9 = set(range(25, 29))
    C_10 = set(range(29, 33))
    P_1 = {9, 10}
    P_2 = {11, 12}
    W_1 = C_N.union(P_1).union(C_9).union(C_10).union(P_2)

    # Based on our analysis, we construct W2, an EJR committee that minimizes s(N,W).
    # This requires 6 candidates from C_N, 10 candidates from C_P (private candidates of voters 1-5),
    # and 4 candidates from A(9) and A(10).
    W_2_C_N = set(range(1, 7))
    W_2_C_P = set(range(9, 19))
    W_2_others = {25, 26, 29, 30}
    W_2 = W_2_C_N.union(W_2_C_P).union(W_2_others)
    
    # Calculate satisfactions
    s_N_W1 = calculate_satisfaction(voter_profiles, group_N, W_1)
    s_N_W2 = calculate_satisfaction(voter_profiles, group_N, W_2)

    # Calculate and print the ratio
    if s_N_W2 == 0:
        ratio = "undefined (division by zero)"
    else:
        ratio = s_N_W1 / s_N_W2
    
    print("s(N,W1) =", s_N_W1)
    print("s(N,W2) =", s_N_W2)
    print(f"The ratio s(N,W1) / s(N,W2) is {s_N_W1} / {s_N_W2} = {ratio}")

solve()
<<<34/29>>>