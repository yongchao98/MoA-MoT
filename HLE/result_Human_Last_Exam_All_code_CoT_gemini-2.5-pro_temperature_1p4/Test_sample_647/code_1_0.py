def solve_committee_problem():
    """
    Solves the committee election problem by calculating the satisfaction
    for two different committees and finding their ratio.
    """

    # Define the approval ballots for the 10 voters.
    # A[i] corresponds to A(i+1)
    A = [
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

    # The group of voters N = {1, ..., 8}. In 0-based indexing, this is {0, ..., 7}.
    N_indices = list(range(8))

    # W1: A core committee with the lowest satisfaction for N.
    # It must contain the 8 common candidates {1..8}.
    # To minimize s(N, W), the other 12 members are chosen from candidates
    # least approved by N: {25..32} (0 approvals) and 4 from {9..24} (1 approval each).
    W1 = set(range(1, 9)) | set(range(25, 33)) | {9, 10, 11, 12}

    # W2: An EJR committee with the lowest satisfaction for N.
    # EJR implies s(N,W) >= 16. This is achieved by selecting the 2 "private"
    # candidates for each voter in N, i.e., {9..24}.
    # The remaining 4 members {25,26,29,30} satisfy EJR for voters 9 and 10.
    W2 = set(range(9, 25)) | {25, 26, 29, 30}

    def calculate_satisfaction(voter_group_indices, committee):
        """Calculates s(N, W) = sum(|A(i) intersect W|) for i in N."""
        total_satisfaction = 0
        for i in voter_group_indices:
            total_satisfaction += len(A[i].intersection(committee))
        return total_satisfaction

    # Calculate the satisfaction for N with W1 and W2.
    s_N_W1 = calculate_satisfaction(N_indices, W1)
    s_N_W2 = calculate_satisfaction(N_indices, W2)

    # Calculate the ratio.
    ratio = s_N_W1 / s_N_W2

    # Print the final equation with all numbers.
    print(f"The satisfaction of group N with committee W1 is s(N,W1) = {s_N_W1}.")
    print(f"The satisfaction of group N with committee W2 is s(N,W2) = {s_N_W2}.")
    print(f"The ratio is {s_N_W1} / {s_N_W2} = {ratio}")

solve_committee_problem()