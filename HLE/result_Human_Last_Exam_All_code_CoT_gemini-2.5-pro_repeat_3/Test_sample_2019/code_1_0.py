def solve_4_queens_toroidal():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard.
    """
    N = 5  # Board size
    M = 4  # Number of queens

    def is_safe(r, c, queens):
        """
        Checks if placing a queen at (r, c) is safe from the existing queens.
        """
        for qr, qc in queens:
            # Check row and column attacks (although row check is redundant in some implementations)
            if r == qr or c == qc:
                return False
            # Check toroidal diagonal attacks
            # Main diagonals: (r - c) mod N is constant
            if (r - c) % N == (qr - qc) % N:
                return False
            # Anti-diagonals: (r + c) mod N is constant
            if (r + c) % N == (qr + qc) % N:
                return False
        return True

    def count_placements(start_index, queens):
        """
        Recursively counts valid queen placements using backtracking.
        - start_index is used to generate combinations instead of permutations, avoiding overcounting.
        - queens is the list of (row, col) coordinates of already placed queens.
        """
        # Base case: If M queens are placed, we found a valid arrangement.
        if len(queens) == M:
            return 1

        count = 0
        # Iterate through all cells on the board, starting from start_index.
        for k in range(start_index, N * N):
            r = k // N
            c = k % N
            
            # Check if it's safe to place a queen here
            if is_safe(r, c, queens):
                # Place the queen
                queens.append((r, c))
                # Recurse to place the next queen, starting from the next cell
                count += count_placements(k + 1, queens)
                # Backtrack: remove the queen to explore other possibilities
                queens.pop()
        
        return count

    # Start the backtracking search from the first cell (index 0) with no queens placed.
    total_ways = count_placements(0, [])
    print(f"The number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard is: {total_ways}")

solve_4_queens_toroidal()
<<<250>>>