import math

def solve_chip_configurations():
    """
    Calculates the number of ways to place 8 chips on an 8x8 board with one chip
    per row and column, such that the placement is symmetric along one of the
    main diagonals.
    """

    # --- Part 1: Configurations symmetric about the main diagonal (|A|) ---
    # This is equivalent to counting the number of involutions on 8 elements (a_8).
    # We use the recurrence relation: a_n = a_{n-1} + (n-1)*a_{n-2}
    n = 8
    a = [0] * (n + 1)
    a[0] = 1
    if n >= 1:
        a[1] = 1
    for i in range(2, n + 1):
        a[i] = a[i - 1] + (i - 1) * a[i - 2]
    n_main_diag = a[n]

    # --- Part 2: Configurations symmetric about the anti-diagonal (|B|) ---
    # This involves pairing 4 row-pairs with 4 column-pairs.
    # The number of ways is 4! * 2^4.
    num_pairs = 4
    n_anti_diag = math.factorial(num_pairs) * (2**num_pairs)

    # --- Part 3: Configurations symmetric about BOTH diagonals (|A intersect B|) ---
    # These permutations (p) must be involutions that commute with s=(1 8)...(4 5).
    # This means p induces an involution (pi) on the 4 pairs of s.
    # We count the number of p based on the structure of pi.
    # Number of involutions on 4 items:
    num_pi_id = 1  # 1 identity pi (4 fixed points, 0 cycles)
    num_pi_one_cycle = math.comb(num_pairs, 2)  # C(4,2)=6 pi with one 2-cycle
    num_pi_two_cycles = math.comb(num_pairs, 2) // 2 # 3 pi with two 2-cycles

    # Contribution from each type of pi:
    # f=fixed points, t=2-cycles. #p = 2^f * 2^t
    # pi=id (f=4, t=0): 1 * 2^4 * 2^0 = 16
    contrib_id = num_pi_id * (2**4 * 2**0)
    # pi has 1 2-cycle (f=2, t=1): 6 * 2^2 * 2^1 = 48
    contrib_one_cycle = num_pi_one_cycle * (2**2 * 2**1)
    # pi has 2 2-cycles (f=0, t=2): 3 * 2^0 * 2^2 = 12
    contrib_two_cycles = num_pi_two_cycles * (2**0 * 2**2)
    
    n_both_diags = contrib_id + contrib_one_cycle + contrib_two_cycles

    # --- Part 4: Final calculation using the Principle of Inclusion-Exclusion ---
    total_configs = n_main_diag + n_anti_diag - n_both_diags
    
    # --- Final Output ---
    print("This problem is solved using the Principle of Inclusion-Exclusion:")
    print("Total = (Symmetric to Main Diagonal) + (Symmetric to Anti-Diagonal) - (Symmetric to Both)")
    print(f"{total_configs} = {n_main_diag} + {n_anti_diag} - {n_both_diags}")

solve_chip_configurations()
<<<1072>>>