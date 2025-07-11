def solve_toroidal_queens():
    """
    Calculates the number of ways to place K non-attacking queens on an N x N toroidal board.
    """
    N = 5  # The size of the board
    K = 4  # The number of queens to place

    def is_safe(r, c, placements):
        """
        Checks if placing a queen at (r, c) is safe from the queens in 'placements'.
        The board is N x N and toroidal.
        """
        for qr, qc in placements:
            # Check row and column attack
            if r == qr or c == qc:
                return False
            # Check toroidal diagonal attacks
            # (a % n + n) % n is used to handle negative results from Python's % operator consistently
            if (r - c - (qr - qc)) % N == 0:
                return False
            if (r + c - (qr + qc)) % N == 0:
                return False
        return True

    def find_placements(queens_to_place, start_cell_index, current_placements):
        """
        A recursive function to find valid placements using backtracking.
        - queens_to_place: Number of queens remaining to be placed.
        - start_cell_index: The board cell (0-24) to start searching from. This
                            ensures we only find combinations, not permutations.
        - current_placements: A list of (row, col) tuples for queens already on the board.
        """
        # Base case: If we have successfully placed all K queens, we've found one solution.
        if queens_to_place == 0:
            return 1

        count = 0
        # Iterate through the board cells starting from 'start_cell_index'.
        for i in range(start_cell_index, N * N):
            # An optimization: if the remaining cells are fewer than queens to place, stop.
            if N * N - i < queens_to_place:
                break
                
            row = i // N
            col = i % N

            # Check if it's safe to place a queen here.
            if is_safe(row, col, current_placements):
                # Place the queen
                current_placements.append((row, col))
                # Recursively call to place the rest of the queens
                count += find_placements(queens_to_place - 1, i + 1, current_placements)
                # Backtrack: remove the queen to explore other possibilities.
                current_placements.pop()
        
        return count

    # Start the search for placing K queens, starting from cell 0 with an empty board.
    total_ways = find_placements(K, 0, [])
    
    print(f"The number of ways to place {K} non-attacking queens on a {N}x{N} toroidal chessboard is:")
    print(total_ways)


solve_toroidal_queens()
<<<550>>>