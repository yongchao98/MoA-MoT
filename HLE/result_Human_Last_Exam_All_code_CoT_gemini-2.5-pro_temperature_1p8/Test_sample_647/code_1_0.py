import math

def calculate_satisfactions_and_ratio():
    """
    Solves the committee election problem by calculating the satisfaction scores
    for the two committees W_1 (core) and W_2 (EJR) and finding their ratio.
    """

    # 1. Define Candidate and Voter Sets
    # A(1)-A(8) approve C_common. Voter i in {1..8} also approves C_unique_i.
    C_common = set(range(1, 9))      # {1, ..., 8}
    C_N_unique = set(range(9, 25))   # {9, ..., 24}
    C_9 = set(range(25, 29))         # {25, ..., 28}
    C_10 = set(range(29, 33))        # {29, ..., 32}
    N = list(range(1, 9)) # Group of voters N = {1, ..., 8}

    # Create approval sets A(i) for voters in N
    A_N = {}
    for i in N:
        # Each voter i in N has two unique candidates: 2*i+7 and 2*i+8
        unique_cands = {2 * i + 7, 2 * i + 8}
        A_N[i] = C_common.union(unique_cands)

    k = 20 # Committee size

    # 2. Calculate s(N, W_1) for the core committee
    # A core committee (maximizing PAV score) will contain all candidates with
    # high approval counts.
    # Analysis shows any core committee W_1 must contain:
    # - All 8 candidates from C_common.
    # - All 8 candidates from C_9 and C_10.
    # This fills 16 slots. The remaining 4 must come from C_N_unique.
    
    num_c_common_in_W1 = 8
    num_c_N_unique_in_W1 = 4
    
    # s(N,W) = sum_{i in N} |A(i) intersect W|
    # For a voter i in N, A(i) = C_common U C_unique_i.
    # A(i) intersect W1 = (C_common intersect W1) U (C_unique_i intersect W1)
    # Since W1 contains all of C_common, |C_common intersect W1| = 8.
    # The satisfaction is the sum of these intersection sizes over all 8 voters in N.
    # Total satisfaction = 8 * (num c_common) + (num c_N_unique)
    s_N_W1 = len(N) * num_c_common_in_W1 + num_c_N_unique_in_W1
    
    # 3. Calculate s(N, W_2) for the EJR committee
    # To find W_2 (EJR committee with lowest satisfaction for N), we must construct
    # an EJR committee that minimizes the number of candidates from C_common.
    # Analysis using the standard EJR definition shows that to satisfy the condition
    # for group N (i.e., there exists a voter i in N with |A(i) intersect W| >= 8),
    # W must contain at least 6 candidates from C_common.
    
    # To minimize s(N, W), we choose the minimum possible number of common candidates.
    num_c_common_in_W2 = 6
    
    # To satisfy EJR for voters 9 and 10, W must have at least 2 members from C_9
    # and 2 from C_10. This accounts for 4 members of W_2.
    num_c_9_10_in_W2 = 4

    # The remaining members of W_2 are 20 - 6 (common) - 4 (from 9,10) = 10.
    # These 10 members must be chosen from C_N_unique.
    num_c_N_unique_in_W2 = k - num_c_common_in_W2 - num_c_9_10_in_W2

    # The satisfaction for N is calculated similarly.
    s_N_W2 = len(N) * num_c_common_in_W2 + num_c_N_unique_in_W2
    
    # 4. Calculate and print the ratio
    ratio = s_N_W1 / s_N_W2
    
    # "Remember in the final code you still need to output each number in the final equation!"
    print(f"The lowest satisfaction for group N with a core committee is s(N,W1) = {s_N_W1}.")
    print(f"The lowest satisfaction for group N with an EJR committee is s(N,W2) = {s_N_W2}.")
    print(f"\nThe ratio is s(N,W1) / s(N,W2) = {s_N_W1} / {s_N_W2} = {ratio}")
    
    # Return the exact fraction and the final numerical value for the answer tag
    frac = math.gcd(s_N_W1, s_N_W2)
    exact_ratio_str = f"{s_N_W1//frac}/{s_N_W2//frac}"
    return ratio, exact_ratio_str

final_ratio, exact_str = calculate_satisfactions_and_ratio()
# Final answer in the required format. Using the exact fraction is best.
# Alternatively, a decimal would be <<<1.1724137931034482>>>
print(f"<<<{exact_str}>>>")
