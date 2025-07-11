def solve_diamond_problem():
    """
    Calculates the largest number of diamonds K such that any valid arrangement
    of K diamonds on an N x N grid guarantees at least one diamond can be moved.
    """
    # The size of the square table
    N = 2024

    # The minimum number of diamonds in an arrangement where no diamond can be moved
    # is M = (N * N) / 4.
    # The question asks for the largest K such that for ANY arrangement of K diamonds,
    # a move is possible. This value is M - 1.
    min_stuck_arrangement_size = (N * N) // 4
    largest_k = min_stuck_arrangement_size - 1

    # Print the final equation and the result
    print(f"The largest value K is calculated by the formula (N * N) / 4 - 1.")
    print(f"For N = {N}:")
    print(f"K = ({N} * {N}) / 4 - 1")
    print(f"K = {N*N} / 4 - 1")
    print(f"K = {min_stuck_arrangement_size} - 1")
    print(f"K = {largest_k}")

solve_diamond_problem()