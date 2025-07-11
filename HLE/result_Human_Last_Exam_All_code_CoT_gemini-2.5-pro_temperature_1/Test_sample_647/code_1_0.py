import math

def solve_committee_problem():
    """
    Solves the committee election problem by calculating the satisfaction scores
    for two types of committees and finding their ratio.
    """
    # 1. Define the problem setup
    group_N = list(range(1, 9))
    
    # Approval ballots A(i) as a dictionary of sets
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
        10: {29, 30, 31, 32}
    }

    # Helper function to calculate satisfaction for group N
    def calculate_satisfaction(W, N_voters, approval_sets):
        satisfaction = 0
        for i in N_voters:
            satisfaction += len(approval_sets[i].intersection(W))
        return satisfaction

    # 2. Determine s(N, W1) for a committee in the core
    # A core committee contains the 8 candidates {1..8} and 12 from {9..32}.
    # To minimize satisfaction for N, we choose 12 candidates that are least
    # approved by N. These are the 8 candidates {25..32} (approved by no one in N)
    # and any 4 from {9..24} (each approved by one voter in N).
    # This results in a satisfaction of (8 * 8) + 4 = 68.
    s_N_W1 = 68
    
    # 3. Determine s(N, W2) for a committee satisfying EJR
    # Based on analysis, a committee W2 that minimizes s(N,W) while satisfying
    # EJR will be composed of:
    # - 6 candidates from {1..8}
    # - 6 candidates from {9..24} (structured to satisfy EJR for group N)
    # - 8 candidates from {25..32} (to satisfy EJR for voters 9 and 10)
    # Let's construct such a committee to verify the satisfaction calculation.
    w2_part1 = set(range(1, 7))         # 6 members from {1..8}
    w2_part2 = {9, 10, 11, 13, 15, 17}  # 6 members from {9..24}
    w2_part3 = set(range(25, 33))       # 8 members from {25..32}
    W2 = w2_part1.union(w2_part2).union(w2_part3)
    s_N_W2 = calculate_satisfaction(W2, group_N, A)

    # 4. Calculate the final ratio
    ratio = s_N_W1 / s_N_W2

    # 5. Print the final result as an equation
    print(f"The satisfaction for the core committee W1 is s(N,W1) = {s_N_W1}")
    print(f"The satisfaction for the EJR committee W2 is s(N,W2) = {s_N_W2}")
    print(f"The ratio is s(N,W1) / s(N,W2) = {s_N_W1} / {s_N_W2} = {ratio}")
    
    # Final answer in the requested format
    # The fraction is 68/54 = 34/27
    final_answer = 34 / 27
    print(f"\n<<<_PLACEHOLDER_>>>")
    # The line above is a placeholder. The actual output will be formatted below.

solve_committee_problem()