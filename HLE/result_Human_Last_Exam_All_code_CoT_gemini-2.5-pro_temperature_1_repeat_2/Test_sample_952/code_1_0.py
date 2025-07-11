def solve_diamond_puzzle():
    """
    Calculates the largest value k for the diamond puzzle.
    """
    # The grid size is 2024x2024
    N = 2024

    # A "stuck" arrangement can be formed by placing diamonds on cells (i, j)
    # where both i and j are even. The number of even indices from 0 to N-1 is N/2.
    side = N // 2

    # The minimum size of a "stuck" arrangement is side * side.
    min_stuck_size = side * side

    # The largest k for which any arrangement is NOT stuck is one less.
    k = min_stuck_size - 1

    # Print the final equation with all its components.
    print(f"The side length for the minimal stuck configuration subgrid is: {side}")
    print(f"The minimum number of diamonds for a 'stuck' arrangement is: {side} * {side} = {min_stuck_size}")
    print(f"The largest value k such that a move is always possible is: {side} * {side} - 1 = {k}")

solve_diamond_puzzle()