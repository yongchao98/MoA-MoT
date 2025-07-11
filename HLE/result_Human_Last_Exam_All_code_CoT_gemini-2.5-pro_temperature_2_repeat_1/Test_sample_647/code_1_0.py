def solve():
    """
    Solves the committee election problem by constructing the committees W1 and W2
    based on analytical derivations and then computing the satisfaction ratio.
    """
    # Define voter approvals as a dictionary of sets
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

    # Group of voters N
    N = set(range(1, 9))

    # Construct W1: core committee with lowest satisfaction for N.
    # This is an AV committee that minimizes satisfaction for N.
    # It includes the 8 most approved candidates {1..8},
    # plus 12 candidates that are worst for N.
    # The 8 candidates {25..32} are not approved by any in N.
    # The next best choice is any 4 from {9..24}. We pick {9,10,11,12}.
    W1 = set(range(1, 9)) | set(range(25, 33)) | {9, 10, 11, 12}

    # Construct W2: EJR committee with lowest satisfaction for N.
    # It must give |A(9) intersect W|>=2, |A(10) intersect W|>=2, and
    # exists i in N with |A(i) intersect W|>=8.
    # To minimize s(N, W), we pick all of {25..32} first.
    # Then we add 12 members from {1..24} to satisfy EJR for N
    # while minimizing satisfaction sum.
    # The minimal choice is 6 from {1..8} and 6 from {9..24}, configured
    # to give one voter satisfaction 8.
    W2 = set(range(25, 33)) | set(range(1, 7)) | {9, 10, 11, 12, 13, 14}


    def calculate_satisfaction(group, committee, approvals):
        """Calculates s(group, committee) = sum_{i in group} |A(i) intersect committee|."""
        total_satisfaction = 0
        for voter_id in group:
            approved_candidates = approvals.get(voter_id, set())
            satisfaction = len(approved_candidates.intersection(committee))
            total_satisfaction += satisfaction
        return total_satisfaction

    s_N_W1 = calculate_satisfaction(N, W1, A)
    s_N_W2 = calculate_satisfaction(N, W2, A)

    ratio = s_N_W1 / s_N_W2

    print(f"W_1 = sorted(list({W1}))")
    print(f"s(N, W_1) = {s_N_W1}")
    print(f"W_2 = sorted(list({W2}))")
    print(f"s(N, W_2) = {s_N_W2}")
    print(f"The ratio s(N,W_1)/s(N,W_2) is {s_N_W1}/{s_N_W2} = {ratio}")

solve()
<<<34/27>>>