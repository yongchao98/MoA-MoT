import math

def solve():
    """
    This function calculates the number of ways to place 8 chips on an 8x8 board
    with one chip per row and column, such that the arrangement is symmetric
    with respect to at least one of the two main diagonals.
    """

    # --- Step 1: Calculate configurations symmetric about the main diagonal (N_main) ---
    # These are permutations that are their own inverse, called involutions.
    # We can count them using the recurrence relation: I_n = I_{n-1} + (n-1)*I_{n-2}.
    # I_0 = 1 (empty set)
    # I_1 = 1 (identity permutation)
    I = [0] * 9
    I[0] = 1
    I[1] = 1
    for n in range(2, 9):
        I[n] = I[n-1] + (n-1) * I[n-2]
    num_main_symm = I[8]

    # --- Step 2: Calculate configurations symmetric about the anti-diagonal (N_anti) ---
    # This corresponds to permutations p where p(9-i) = 9-p(i).
    # The problem breaks down into pairing rows {1..4} with columns.
    # There are 4 pairs of rows {i, 9-i} and 4 pairs of columns {j, 9-j}.
    # We permute the column pairs (4! ways) and for each pair, there are 2 ways to assign columns.
    n_pairs = 8 // 2
    num_anti_symm = math.factorial(n_pairs) * (2**n_pairs)

    # --- Step 3: Calculate configurations symmetric about BOTH diagonals (N_both) ---
    # These are involutions that also satisfy the anti-diagonal symmetry condition.
    # This calculation depends on the cycle structure of involutions on the 4 pairs of indices.
    # An involution q on 4 items can have:
    # 1. Four 1-cycles (identity): 1 case. Contributes 1 * 2^4 configurations.
    # 2. One 2-cycle and two 1-cycles: C(4,2) = 6 cases. Each contributes 6 * 2^3 configurations.
    # 3. Two 2-cycles: 3 cases (C(4,2)*C(2,2)/2). Each contributes 3 * 2^2 configurations.
    
    # Case 1: 0 two-cycles, 4 one-cycles (q = id)
    # Number of such involutions q on 4 pairs
    num_q_0_twocycles = 1 
    # Number of p for each q
    num_p_for_q0 = 2**4 
    
    # Case 2: 1 two-cycle, 2 one-cycles
    num_q_1_twocycle = math.comb(4, 2)
    num_p_for_q1 = 2**3 # Cycles in q: (1,2)(3)(4), 3 cycles

    # Case 3: 2 two-cycles
    num_q_2_twocycles = math.comb(4, 2) * math.comb(2, 2) // 2
    num_p_for_q2 = 2**2 # Cycles in q: (1,2)(3,4), 2 cycles
    
    num_both_symm = (num_q_0_twocycles * num_p_for_q0) + \
                    (num_q_1_twocycle * num_p_for_q1) + \
                    (num_q_2_twocycles * num_p_for_q2)

    # --- Step 4: Apply Principle of Inclusion-Exclusion ---
    total_configs = num_main_symm + num_anti_symm - num_both_symm
    
    print(f"Number of configurations symmetric about the main diagonal: {num_main_symm}")
    print(f"Number of configurations symmetric about the anti-diagonal: {num_anti_symm}")
    print(f"Number of configurations symmetric about both diagonals: {num_both_symm}")
    print("\nApplying the Principle of Inclusion-Exclusion:")
    print(f"Total possible configurations = (Main Diagonal) + (Anti-Diagonal) - (Both)")
    print(f"Total = {num_main_symm} + {num_anti_symm} - {num_both_symm} = {total_configs}")

solve()