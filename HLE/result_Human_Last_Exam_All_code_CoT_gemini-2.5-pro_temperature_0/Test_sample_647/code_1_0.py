import math

def calculate_satisfaction(group_N, committee_W, all_ballots):
    """Calculates the satisfaction of a group N with a committee W."""
    total_satisfaction = 0
    for voter_id in group_N:
        # Voter ballots are 1-indexed, Python lists are 0-indexed
        voter_ballot = all_ballots[voter_id - 1]
        satisfaction_i = len(voter_ballot.intersection(committee_W))
        total_satisfaction += satisfaction_i
    return total_satisfaction

def solve():
    """
    Solves the committee election problem.
    """
    # Define the approval ballots for all 10 voters
    A = [
        {1,2,3,4,5,6,7,8,9,10},
        {1,2,3,4,5,6,7,8,11,12},
        {1,2,3,4,5,6,7,8,13,14},
        {1,2,3,4,5,6,7,8,15,16},
        {1,2,3,4,5,6,7,8,17,18},
        {1,2,3,4,5,6,7,8,19,20},
        {1,2,3,4,5,6,7,8,21,22},
        {1,2,3,4,5,6,7,8,23,24},
        {25,26,27,28},
        {29,30,31,32}
    ]

    # Define the group of voters N
    N = range(1, 9) # Voters 1 to 8

    # --- Step 1: Determine W1 (Core committee) and s(N, W1) ---
    # A committee in the core is a PAV-optimal committee.
    # Analysis shows that to maximize the PAV score, we must select:
    # 1. All 8 candidates shared by voters in N: {1..8}
    # 2. All 8 candidates approved by voters 9 and 10: {25..32}
    # 3. The remaining 4 candidates (to reach size 20) must come from the set
    #    of candidates uniquely approved by one voter in N: {9..24}.
    # To find the lowest satisfaction for N, any choice of these 4 is equivalent,
    # as each adds exactly 1 to the satisfaction of one voter in N.
    
    C_shared = set(range(1, 9))
    C_9_10 = set(range(25, 33))
    # We pick 4 candidates from C_unique_N, e.g., {9, 10, 11, 12}
    C_unique_N_sample = set(range(9, 13))
    
    W1 = C_shared.union(C_9_10).union(C_unique_N_sample)
    
    s_N_W1 = calculate_satisfaction(N, W1, A)
    
    print("--- Analysis for W1 (Core) ---")
    print(f"The committee W1 is composed of:")
    print(f" - 8 candidates shared by group N: {sorted(list(C_shared))}")
    print(f" - 8 candidates from voters 9 and 10: {sorted(list(C_9_10))}")
    print(f" - 4 candidates from the remaining pool: {sorted(list(C_unique_N_sample))}")
    print(f"The total satisfaction for group N with W1 is s(N,W1) = {s_N_W1}")
    print("-" * 20)

    # --- Step 2: Determine W2 (EJR committee) and s(N, W2) ---
    # To satisfy EJR while minimizing s(N,W), we construct W2 as follows:
    # 1. Maximize candidates from outside N's approval sets: take all 8 from {25..32}.
    # 2. This leaves 12 spots to fill from {1..24}.
    # 3. EJR for group N requires at least one voter i in N to have |A(i) intersect W2| >= 8.
    # 4. To achieve this with minimum s(N,W), we should pick 6 candidates from C_shared {1..8}
    #    and 2 unique candidates for one voter (e.g., {9,10} for voter 1).
    #    This gives voter 1 a satisfaction of 6+2=8.
    # 5. The remaining 12 - (6+2) = 4 candidates are chosen from the rest of {9..24}
    #    to complete the committee.
    
    W2_C_9_10 = set(range(25, 33)) # 8 candidates
    W2_C_shared = set(range(1, 7)) # 6 candidates
    # 6 candidates from C_unique_N to satisfy EJR for N and fill the committee
    W2_C_unique_N = set(range(9, 15)) # {9,10} for voter 1, {11,12} for v2, {13,14} for v3
    
    W2 = W2_C_9_10.union(W2_C_shared).union(W2_C_unique_N)

    s_N_W2 = calculate_satisfaction(N, W2, A)

    print("--- Analysis for W2 (EJR) ---")
    print(f"The committee W2 is composed of:")
    print(f" - 8 candidates from voters 9 and 10: {sorted(list(W2_C_9_10))}")
    print(f" - 6 candidates shared by group N: {sorted(list(W2_C_shared))}")
    print(f" - 6 candidates from the remaining pool: {sorted(list(W2_C_unique_N))}")
    print(f"The total satisfaction for group N with W2 is s(N,W2) = {s_N_W2}")
    print("-" * 20)

    # --- Step 3: Calculate the ratio ---
    ratio = s_N_W1 / s_N_W2
    
    print("--- Final Calculation ---")
    print(f"The ratio s(N,W1) / s(N,W2) is:")
    print(f"{s_N_W1} / {s_N_W2} = {ratio}")
    
    return ratio

# Run the solver and print the final answer in the required format
final_ratio = solve()
print(f"\n<<<{s_N_W1 / s_N_W2}>>>")