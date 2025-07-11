def solve_toroidal_queens():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard.
    """
    BOARD_SIZE = 5
    NUM_QUEENS = 4

    def is_safe(r2, c2, placed_queens):
        """
        Checks if placing a new queen at (r2, c2) is safe with respect to the already placed queens.
        A placement is safe if the new queen is not on the same row, column, or toroidal diagonal
        as any existing queen.
        """
        for r1, c1 in placed_queens:
            # Check for row or column attack
            if r1 == r2 or c1 == c2:
                return False
            # Check for toroidal diagonal attacks using modular arithmetic
            # Main diagonals: (r1 - c1) mod N == (r2 - c2) mod N
            if (r1 - c1) % BOARD_SIZE == (r2 - c2) % BOARD_SIZE:
                return False
            # Anti-diagonals: (r1 + c1) mod N == (r2 + c2) mod N
            if (r1 + c1) % BOARD_SIZE == (r2 + c2) % BOARD_SIZE:
                return False
        return True

    def count_placements(k, start_index, placed_queens):
        """
        Recursively counts valid queen placements using backtracking.

        :param k: The number of queens remaining to be placed.
        :param start_index: The linear index on the board to start searching from. This
                            ensures we only count combinations, not permutations.
        :param placed_queens: A list of (row, col) tuples for queens already placed.
        :return: The number of valid placements found from this state.
        """
        # Base case: If k is 0, we have successfully placed all queens.
        if k == 0:
            return 1

        count = 0
        # Iterate through all possible squares starting from start_index.
        for i in range(start_index, BOARD_SIZE * BOARD_SIZE):
            r = i // BOARD_SIZE
            c = i % BOARD_SIZE

            # Check if placing a queen at (r, c) is safe.
            if is_safe(r, c, placed_queens):
                # If safe, place the queen and recurse.
                placed_queens.append((r, c))
                count += count_placements(k - 1, i + 1, placed_queens)
                # Backtrack: remove the queen to explore other possibilities.
                placed_queens.pop()
        
        return count

    # Start the calculation.
    total_ways = count_placements(NUM_QUEENS, 0, [])
    print(total_ways)

solve_toroidal_queens()
<<<200>>>