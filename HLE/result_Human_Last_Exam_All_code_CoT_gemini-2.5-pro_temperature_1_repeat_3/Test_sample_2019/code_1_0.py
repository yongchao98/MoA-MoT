def solve_queens_toroidal():
    """
    Calculates the number of ways to place a given number of non-attacking queens
    on a toroidal chessboard of a given size.
    """
    BOARD_SIZE = 5
    NUM_QUEENS = 4

    def solve(k, start_index, queens):
        """
        Recursively finds the number of ways to place k non-attacking queens.

        Args:
            k: The number of queens remaining to be placed.
            start_index: The linear index on the board to start searching from. This ensures
                         that we count combinations of queen placements, not permutations.
            queens: A list of (row, col) tuples for queens already placed.

        Returns:
            The number of valid placements found from this state.
        """
        # Base case: If we have successfully placed all queens, we've found one solution.
        if k == 0:
            return 1

        count = 0
        # Iterate through possible squares for the next queen.
        # The upper bound is optimized: we need to leave at least (k-1) squares
        # for the remaining queens.
        for i in range(start_index, BOARD_SIZE * BOARD_SIZE - (k - 1)):
            r = i // BOARD_SIZE
            c = i % BOARD_SIZE

            # Check if placing a queen at (r, c) is safe with respect to existing queens.
            is_placement_safe = True
            for qr, qc in queens:
                # Check for row or column attacks.
                if r == qr or c == qc:
                    is_placement_safe = False
                    break
                
                # Check for toroidal diagonal attacks.
                dr = abs(r - qr)
                dc = abs(c - qc)
                
                toroidal_dr = min(dr, BOARD_SIZE - dr)
                toroidal_dc = min(dc, BOARD_SIZE - dc)
                
                if toroidal_dr == toroidal_dc:
                    is_placement_safe = False
                    break
            
            if is_placement_safe:
                # If safe, place the queen and recurse.
                queens.append((r, c))
                count += solve(k - 1, i + 1, queens)
                # Backtrack: remove the queen to explore other possibilities.
                queens.pop()
                
        return count

    # Start the backtracking search.
    total_ways = solve(NUM_QUEENS, 0, [])
    
    print(f"The number of ways to place {NUM_QUEENS} non-attacking queens on a {BOARD_SIZE}x{BOARD_SIZE} toroidal chessboard is: {total_ways}")

solve_queens_toroidal()