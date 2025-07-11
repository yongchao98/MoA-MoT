def calculate_satisfaction(approvals, committee):
    """Calculates the satisfaction of a group of voters with a committee."""
    total_satisfaction = 0
    for voter_approvals in approvals:
        # |A(i) intersect W| is the size of the intersection of two sets
        intersection = voter_approvals.intersection(committee)
        total_satisfaction += len(intersection)
    return total_satisfaction

def solve():
    """
    Solves the committee election problem by determining the committees W1 and W2,
    calculating their respective satisfactions for group N, and finding the ratio.
    """
    # Define the approval sets for all 10 voters
    all_approvals = [
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

    # Group N consists of the first 8 voters
    group_n_approvals = all_approvals[0:8]

    # --- Step 1: Determine W1 and s(N, W1) ---
    # W1 is a committee in the core with the lowest satisfaction for N.
    # Our analysis suggests such a committee will give 16 proportionally-deserved seats to group N.
    # To be stable against blocking by N, it must contain all 8 of N's common candidates {1..8}.
    # The other 8 seats for N come from their "private" candidates {9..24}.
    # The remaining 4 seats go to voters 9 and 10 (2 each).
    # Any such committee will yield the same satisfaction for N.
    
    # Example W1 construction:
    # 8 common candidates for N
    w1_common_n = set(range(1, 9))
    # 8 private candidates for N (e.g., one from each voter's private pair)
    w1_private_n = set(range(9, 17))
    # 2 candidates for voter 9
    w1_v9 = {25, 26}
    # 2 candidates for voter 10
    w1_v10 = {29, 30}
    W1 = w1_common_n.union(w1_private_n).union(w1_v9).union(w1_v10)
    
    s_n_w1 = calculate_satisfaction(group_n_approvals, W1)

    # --- Step 2: Determine W2 and s(N, W2) ---
    # W2 is an EJR committee with the lowest satisfaction for N.
    # EJR requires that for group N, some voter j must have satisfaction >= 8.
    # To minimize total satisfaction for N, we build W2 to satisfy this constraint minimally.
    # We give 6 common candidates and 2 private candidates (for one voter, e.g., voter 1).
    # The remaining seats are filled to avoid adding more satisfaction to N.
    
    # W2 construction:
    # 6 common candidates from {1..8}
    w2_common_n = set(range(1, 7))
    # 2 private candidates to satisfy EJR for voter 1 (making their sat 6+2=8)
    w2_ejr_sats = {9, 10}
    # To minimize N's satisfaction, fill remaining seats from outside N's approval set union
    w2_others = set(range(25, 33)) # 8 seats for voters 9 and 10
    # Fill the final 4 seats from other private candidates of N
    w2_fillers = set(range(11, 15))
    W2 = w2_common_n.union(w2_ejr_sats).union(w2_others).union(w2_fillers)
    
    s_n_w2 = calculate_satisfaction(group_n_approvals, W2)

    # --- Step 3: Calculate and print the ratio ---
    ratio = s_n_w1 / s_n_w2

    print(f"The satisfaction of N with W1, s(N,W1), is {s_n_w1}.")
    print(f"The satisfaction of N with W2, s(N,W2), is {s_n_w2}.")
    print("The required ratio is s(N,W1) / s(N,W2):")
    # The prompt requires printing each number in the final equation.
    print(f"{s_n_w1} / {s_n_w2} = {ratio}")
    return ratio

# Execute the solution
final_ratio = solve()
# Present the final numerical answer in the required format
print(f"\n<<<{final_ratio}>>>")