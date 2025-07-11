def solve_toroidal_queens():
    """
    Calculates the number of ways to place K non-attacking queens on an N x N toroidal chessboard.
    """
    N = 5  # Size of the chessboard
    K = 4  # Number of queens to place

    def is_safe(r2, c2, queens):
        """
        Checks if placing a queen at (r2, c2) is safe from the list of existing queens.
        'queens' is a list of (row, col) tuples.
        """
        for r1, c1 in queens:
            # Check for row or column attack (should not happen with the main loop's logic, but good for safety)
            if r1 == r2 or c1 == c2:
                return False

            # Check for toroidal diagonal attack
            row_dist = min(abs(r1 - r2), N - abs(r1 - r2))
            col_dist = min(abs(c1 - c2), N - abs(c1 - c2))
            if row_dist == col_dist:
                return False
        return True

    def count_placements(k, start_cell, placed_queens):
        """
        Recursively counts the number of valid placements.
        - k: number of queens left to place
        - start_cell: the starting cell index (0-24) for the search to avoid duplicates
        - placed_queens: a list of (r, c) tuples for queens already on the board
        """
        # Base case: if all K queens are placed, we found one valid solution.
        if k == 0:
            return 1

        total_count = 0
        # Iterate through all possible cells starting from start_cell.
        # This ensures we only generate combinations of queen placements, not permutations.
        for i in range(start_cell, N * N):
            r = i // N
            c = i % N

            # Check if placing a queen at (r, c) is safe
            if is_safe(r, c, placed_queens):
                # If it's safe, place the queen and recurse to place the next one.
                # The search for the next queen starts from the next cell (i + 1).
                new_queens = placed_queens + [(r, c)]
                total_count += count_placements(k - 1, i + 1, new_queens)
        
        return total_count

    # Start the calculation
    total_ways = count_placements(K, 0, [])

    print(f"Board size: {N}x{N}")
    print(f"Number of queens: {K}")
    print(f"The number of ways to place {K} non-attacking queens on a {N}x{N} toroidal chessboard is: {total_ways}")
    return total_ways

result = solve_toroidal_queens()
# The final answer format is required by the problem statement.
# <<<answer>>>
# where answer is the number computed.
final_answer_str = f"<<<{result}>>>"