def solve_diamond_problem():
    """
    Calculates the largest value K based on the logic derived.

    The problem is to find the largest integer K such that any valid arrangement
    of K diamonds on a 2024x2024 grid has at least one "movable" diamond.

    This is equivalent to finding M-1, where M is the minimum number of
    diamonds required to form a "stuck" arrangement (where no diamond is movable).

    We can construct a stuck arrangement by placing diamonds on cells (r, c)
    where both r and c are even (using 0-based indexing).
    For a 2024x2024 grid, the even indices are 0, 2, ..., 2022.

    The number of choices for an even row index is 2024 / 2 = 1012.
    The number of choices for an even column index is 2024 / 2 = 1012.

    So, the size of this stuck arrangement is M = 1012 * 1012.
    This is known to be the minimum size for a stuck arrangement.

    The value K is therefore M - 1.
    """
    grid_size = 2024

    # Number of even indices from 0 to grid_size - 1
    half_size = grid_size // 2

    # M is the minimum number of diamonds in a "stuck" configuration
    M = half_size * half_size

    # K is the largest number of diamonds for which a move is always possible
    K = M - 1

    # The problem asks to output each number in the final equation.
    print(f"The grid size is {grid_size}x{grid_size}.")
    print(f"The minimum number of diamonds in a 'stuck' configuration (M) is found by constructing a pattern.")
    print(f"The pattern places diamonds on all cells (row, col) where row and col are even.")
    print(f"Number of even rows/columns = {grid_size} / 2 = {half_size}.")
    print(f"M = {half_size} * {half_size} = {M}")
    print(f"The largest value K such that any arrangement of K diamonds has a movable piece is M - 1.")
    print(f"K = {M} - 1 = {K}")
    print("\nFinal equation:")
    print(f"{half_size} * {half_size} - 1 = {K}")


solve_diamond_problem()