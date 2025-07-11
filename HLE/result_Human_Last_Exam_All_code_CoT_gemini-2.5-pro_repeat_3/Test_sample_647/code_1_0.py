def solve_committee_problem():
    """
    Calculates the satisfaction values for committees W1 (core) and W2 (EJR)
    and computes their ratio.
    """

    # --- Step 1: Calculate s(N, W_1), the lowest satisfaction for group N in a core committee ---

    # In a core committee, seats are allocated proportionally to group size.
    # Group N (8 voters) quota: (8/10) * 20 = 16 seats.
    # To minimize satisfaction s(N, W) = 8 * (num_common) + 1 * (num_unique),
    # we should pick 0 common candidates and 16 unique candidates for group N.
    num_common_w1 = 0
    num_unique_n_w1 = 16
    s_N_W1 = 8 * num_common_w1 + 1 * num_unique_n_w1
    print(f"For the core committee W_1 with the lowest satisfaction for N:")
    print(f"  Number of common candidates chosen = {num_common_w1}")
    print(f"  Number of unique candidates for N chosen = {num_unique_n_w1}")
    print(f"s(N,W_1) = 8 * {num_common_w1} + 1 * {num_unique_n_w1} = {s_N_W1}")
    print("-" * 20)

    # --- Step 2: Calculate s(N, W_2), the lowest satisfaction for group N in an EJR committee ---

    # To minimize s(N, W) for an EJR committee, we select as many candidates
    # as possible that are NOT approved by group N, while still meeting EJR guarantees.
    # We can select all 8 candidates approved by voters 9 and 10.
    # This leaves 20 - 8 = 12 seats to fill.
    # To keep s(N, W) low, we fill these 12 seats with unique candidates for N.
    num_common_w2 = 0
    num_unique_n_w2 = 12
    s_N_W2 = 8 * num_common_w2 + 1 * num_unique_n_w2
    print(f"For the EJR committee W_2 with the lowest satisfaction for N:")
    print(f"  Number of common candidates chosen = {num_common_w2}")
    print(f"  Number of unique candidates for N chosen = {num_unique_n_w2}")
    print(f"s(N,W_2) = 8 * {num_common_w2} + 1 * {num_unique_n_w2} = {s_N_W2}")
    print("-" * 20)

    # --- Step 3: Calculate and print the final ratio ---
    ratio = s_N_W1 / s_N_W2
    print("The final ratio is:")
    print(f"s(N,W_1) / s(N,W_2) = {s_N_W1} / {s_N_W2} = {ratio}")

    # Final answer in the required format
    print(f"\n<<<{ratio}>>>")

solve_committee_problem()