def solve_diamond_problem():
    """
    Calculates the largest value k for the diamond arrangement problem.

    The problem asks for the largest integer k such that for any valid arrangement of k
    diamonds on a 2024x2024 grid, there is always one diamond that can be moved to an
    adjacent cell while keeping the arrangement valid.

    This is equivalent to finding M-1, where M is the size of the smallest "stuck"
    arrangement (an independent 2-dominating set).

    For an n x n grid, the size M is given by a known formula based on n mod 6.
    - M = floor(n^2 / 3)      if n mod 6 in {0, 1, 3, 4}
    - M = floor(n^2 / 3) + 1  if n mod 6 in {2, 5}
    """
    n = 2024

    # Determine which formula to use based on n mod 6
    n_mod_6 = n % 6

    # Calculate n^2
    n_squared = n * n

    # Calculate floor(n^2 / 3)
    floor_val = n_squared // 3

    # Calculate M, the size of the smallest stuck arrangement
    if n_mod_6 in [2, 5]:
        M = floor_val + 1
    else:
        M = floor_val
    
    # The answer k is M - 1
    k = M - 1

    print(f"Grid size n = {n}")
    print(f"n^2 = {n_squared}")
    print(f"n mod 6 = {n_mod_6}, so we use the formula M = floor(n^2 / 3) + 1")
    print(f"M = floor({n_squared} / 3) + 1 = {floor_val} + 1 = {M}")
    print(f"The largest value k is M - 1.")
    print(f"k = {M} - 1 = {k}")


solve_diamond_problem()