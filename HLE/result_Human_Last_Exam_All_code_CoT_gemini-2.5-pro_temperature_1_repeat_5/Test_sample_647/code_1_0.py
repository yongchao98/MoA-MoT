import math

def solve_committee_election():
    """
    This function solves the committee election problem by constructing the committees
    W_1 and W_2 and calculating the ratio of their satisfaction for group N.
    """
    
    # Define the approval sets for all 10 voters
    A = {
        1: set(range(1, 11)),
        2: set(range(1, 9)) | {11, 12},
        3: set(range(1, 9)) | {13, 14},
        4: set(range(1, 9)) | {15, 16},
        5: set(range(1, 9)) | {17, 18},
        6: set(range(1, 9)) | {19, 20},
        7: set(range(1, 9)) | {21, 22},
        8: set(range(1, 9)) | {23, 24},
        9: set(range(25, 29)),
        10: set(range(29, 33)),
    }
    
    # Define the group of voters N
    N = set(range(1, 9))
    
    # --- Determine W_1 and calculate s(N, W_1) ---
    # W_1 is a core committee that minimizes satisfaction for N.
    # Utilitarian committees are in the core. They contain the 8 candidates
    # {1..8} approved by most voters.
    # To minimize s(N,W), the other 12 members are chosen to be as
    # disagreeable to voters in N as possible.
    # This means choosing all 8 candidates from voters 9 and 10 ({25..32})
    # and the remaining 4 from the set {9..24}.
    
    W1_common_core = set(range(1, 9))
    W1_from_others = set(range(25, 33))
    W1_from_N_unique = set(range(9, 13)) # Example: 4 candidates from {9..24}
    
    W1 = W1_common_core | W1_from_others | W1_from_N_unique
    
    s_N_W1 = 0
    for i in N:
        s_N_W1 += len(A[i].intersection(W1))

    # --- Determine W_2 and calculate s(N, W_2) ---
    # W_2 is an EJR committee that minimizes satisfaction for N.
    # Analysis of the EJR constraints shows this is achieved by a committee that ensures
    # at least 4 voters in N have satisfaction >= 8, while keeping the overall
    # count of candidates approved by N as low as possible.
    # This committee contains:
    # - 6 candidates from {1..8}
    # - The 8 unique candidates of 4 voters from N (to raise their satisfaction to 8)
    # - The remaining 6 candidates from those approved by voters 9 and 10
    
    W2_common_part = set(range(1, 7))
    W2_N_unique_part = set(range(9, 17)) # Candidates from B_1, B_2, B_3, B_4
    W2_from_others = set(range(25, 31))  # 6 candidates from {25..32}
    
    W2 = W2_common_part | W2_N_unique_part | W2_from_others
    
    s_N_W2 = 0
    for i in N:
        s_N_W2 += len(A[i].intersection(W2))

    # --- Calculate and print the final ratio ---
    ratio = s_N_W1 / s_N_W2
    
    print(f"Satisfaction for group N with committee W_1: s(N,W_1) = {s_N_W1}")
    print(f"Satisfaction for group N with committee W_2: s(N,W_2) = {s_N_W2}")
    print(f"The ratio s(N,W_1) / s(N,W_2) is {s_N_W1} / {s_N_W2} = {ratio}")

solve_committee_election()
print(f'<<<{68/56}>>>')