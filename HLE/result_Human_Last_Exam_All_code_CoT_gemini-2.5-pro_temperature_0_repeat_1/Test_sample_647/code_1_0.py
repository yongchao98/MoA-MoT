def solve_committee_problem():
    """
    Solves the committee election problem by analyzing EJR and Core properties.
    """
    # --- Problem Parameters ---
    num_voters = 10
    committee_size_k = 20
    # Candidate sets are defined implicitly by the logic below.
    # C_N_common: {1..8}, C_N_unique: {9..24}, C_9: {25..28}, C_10: {29..32}
    
    print("Analyzing the problem to find the satisfactions for W1 and W2.")
    print("-" * 50)

    # --- Step 1 & 2: Analysis of EJR and finding s(N, W_2) ---
    print("Analysis for W_2 (lowest satisfaction for N under EJR):")
    
    # The EJR rule states that for any group of voters N' of size |N'| >= j * (n/k)
    # that is j-cohesive (all approve at least j common candidates),
    # they must be represented by at least j of those candidates in the committee.
    # Here, n=10, k=20, so n/k = 0.5.

    # For group N = {1..8}, |N|=8. Let j=8.
    # Condition 1: |N| >= j * n/k  =>  8 >= 8 * 0.5  =>  8 >= 4. (Holds)
    # Condition 2: Group N is 8-cohesive as they all approve {1..8}. (Holds)
    # EJR implies the committee must contain all 8 of these common candidates.
    k_common_in_W2 = 8
    
    # For group {9}, |{9}|=1. Let j=2.
    # Condition 1: |{9}| >= j * n/k => 1 >= 2 * 0.5 => 1 >= 1. (Holds)
    # Condition 2: Group {9} is 4-cohesive, which is >= 2. (Holds)
    # EJR implies the committee must contain at least 2 candidates from A(9).
    
    # Similarly for group {10}, the committee must contain at least 2 from A(10).

    # The satisfaction for group N is s(N, W) = 8 * |W intersect C_N_common| + |W intersect C_N_unique|.
    # Since an EJR committee W must contain C_N_common, this becomes s(N, W) = 8*8 + |W intersect C_N_unique|.
    # s(N, W) = 64 + k_unique.
    
    # To minimize satisfaction, we must minimize k_unique.
    # The committee size is 20: k_common + k_unique + k_9 + k_10 = 20.
    # With k_common=8, we have: k_unique + k_9 + k_10 = 12.
    # To minimize k_unique, we must maximize k_9 and k_10, subject to EJR constraints (k_9>=2, k_10>=2).
    # We can choose the maximum possible values: k_9=4 and k_10=4.
    k_unique_in_W2 = 12 - 4 - 4
    
    s_N_W2 = 64 + k_unique_in_W2
    print(f"An EJR committee W_2 minimizing s(N,W) must have {k_unique_in_W2} candidates from N's unique approvals.")
    print(f"The satisfaction s(N, W_2) = 8 * {k_common_in_W2} + {k_unique_in_W2} = {s_N_W2}")
    print("-" * 50)

    # --- Step 3 & 4: Analysis of the Core and finding s(N, W_1) ---
    print("Analysis for W_1 (lowest satisfaction for N in the Core):")

    # A committee is in the core if no group of voters can block it.
    # Group N={1..8} is very cohesive. If a committee W does not contain all 8 of their
    # common candidates, N can propose an alternative committee W' that *does* contain
    # all 8. This makes every voter in N strictly better off, blocking W.
    # Therefore, any committee in the core must contain all 8 common candidates of N.
    k_common_in_W1 = 8

    # The logic to find the minimum satisfaction is now identical to the EJR case.
    # s(N, W) = 64 + k_unique.
    # To minimize k_unique, we must maximize the number of candidates from outside N's approval sets.
    # k_unique + k_9 + k_10 = 12.
    # The core property does not impose minimums on k_9 or k_10, so we can maximize them to their limit.
    k_unique_in_W1 = 12 - 4 - 4

    s_N_W1 = 64 + k_unique_in_W1
    print(f"A core committee W_1 minimizing s(N,W) must have {k_unique_in_W1} candidates from N's unique approvals.")
    print(f"The satisfaction s(N, W_1) = 8 * {k_common_in_W1} + {k_unique_in_W1} = {s_N_W1}")
    print("-" * 50)

    # --- Step 5: Final Calculation ---
    print("Final Calculation:")
    if s_N_W2 == 0:
        print("s(N, W_2) is zero, division is not possible.")
        ratio = "undefined"
    else:
        ratio = s_N_W1 / s_N_W2
    
    print(f"The ratio s(N,W_1)/s(N,W_2) is calculated as: {s_N_W1} / {s_N_W2} = {ratio}")
    print(f"<<<{ratio}>>>")

solve_committee_problem()