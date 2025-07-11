def solve_diamond_problem():
    """
    Calculates the largest value K such that for any valid arrangement of K diamonds,
    at least one can be moved to an adjacent cell.
    """
    # The size of the square table
    N = 2024

    # The problem is equivalent to finding M-1, where M is the minimum number
    # of diamonds in a "stuck" configuration.
    # For an N x N grid where N is even, M = (N/2) * (N/2).
    half_N = N // 2

    # M is the minimum number of diamonds for a stuck arrangement
    M = half_N * half_N

    # K is the largest number of diamonds for which no arrangement is stuck
    K = M - 1

    # Print the final equation with all numbers, as requested.
    print(f"The minimum number of diamonds for a 'stuck' arrangement (M) on a {N}x{N} grid is:")
    print(f"M = ({N} / 2) * ({N} / 2) = {half_N} * {half_N} = {M}")
    print("\nThe largest value K such that any arrangement is NOT stuck is M - 1:")
    print(f"K = {M} - 1 = {K}")

solve_diamond_problem()