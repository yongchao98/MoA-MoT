import math

def solve_chip_configurations():
    """
    This function calculates the number of possible symmetric chip configurations on an 8x8 board.
    """
    
    # --- Step 1: Calculate N_main ---
    # N_main is the number of involutions on 8 elements, a_8.
    # Using the recurrence a_n = a_{n-1} + (n-1) * a_{n-2}
    # a_0 = 1, a_1 = 1
    a = [0] * 9
    a[0] = 1
    a[1] = 1
    for n in range(2, 9):
        a[n] = a[n-1] + (n-1) * a[n-2]
    n_main = a[8]

    # --- Step 2: Calculate N_anti ---
    # For anti-diagonal symmetry, placing a chip in row i determines the chip in row 9-i.
    # Choices for p(1): 8, which determines p(8).
    # Choices for p(2): 6, which determines p(7).
    # Choices for p(3): 4, which determines p(6).
    # Choices for p(4): 2, which determines p(5).
    n_anti = 8 * 6 * 4 * 2

    # --- Step 3: Calculate N_both ---
    # For symmetry about both diagonals, the permutation p must be an involution
    # and satisfy p(9-i) = 9-p(i).
    # This structure partitions the indices {1..8} into 4 pairs S_i = {i, 9-i}.
    # The permutation p induces an involution on these 4 pairs.
    
    # Case 1: 0 swaps (identity involution on pairs).
    # Number of ways: 1 * 2^4 = 16
    n_both_c1 = 2**4

    # Case 2: 1 swap (e.g., S1 <-> S2, S3, S4 fixed).
    # Number of such involutions on pairs: C(4,2) = 6.
    # For each, there are 2 (for the swap) * 2^2 (for fixed) = 8 ways.
    # Number of ways: 6 * 8 = 48
    num_pair_involutions_1_swap = math.comb(4, 2)
    n_both_c2 = num_pair_involutions_1_swap * (2 * 2**2)

    # Case 3: 2 swaps (e.g., S1 <-> S2, S3 <-> S4).
    # Number of such involutions on pairs: C(4,2)/2 = 3.
    # For each, there are 2 * 2 = 4 ways.
    # Number of ways: 3 * 4 = 12
    num_pair_involutions_2_swaps = math.comb(4, 2) // 2
    n_both_c3 = num_pair_involutions_2_swaps * (2 * 2)
    
    n_both = n_both_c1 + n_both_c2 + n_both_c3

    # --- Step 4: Final Calculation ---
    total_configs = n_main + n_anti - n_both

    print(f"Number of configurations symmetric along the main diagonal (N_main): {n_main}")
    print(f"Number of configurations symmetric along the anti-diagonal (N_anti): {n_anti}")
    print(f"Number of configurations symmetric along both diagonals (N_both): {n_both}")
    print("\nUsing the Principle of Inclusion-Exclusion, the total number of configurations is:")
    print(f"Total = N_main + N_anti - N_both")
    print(f"Total = {n_main} + {n_anti} - {n_both} = {total_configs}")

solve_chip_configurations()
print("\n<<<1072>>>")