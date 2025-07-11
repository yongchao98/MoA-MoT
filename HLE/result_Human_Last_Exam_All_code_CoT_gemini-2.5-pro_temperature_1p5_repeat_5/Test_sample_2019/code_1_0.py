def solve_4_queens_toroidal():
    """
    Calculates the number of ways to place 4 non-attacking queens on a 5x5 toroidal chessboard.
    """
    board_size = 5
    num_queens = 4
    count = 0

    def backtrack(k, start_index, positions):
        """
        Recursively places queens on the board.

        :param k: The number of queens placed so far.
        :param start_index: The starting cell index for placing the next queen.
        :param positions: A list of (row, col) tuples for the queens already placed.
        """
        nonlocal count
        # Base case: if all queens are placed, we found a valid solution.
        if k == num_queens:
            count += 1
            return

        # Pruning: if remaining cells are not enough to place remaining queens, stop.
        if board_size * board_size - start_index < num_queens - k:
            return

        # Iterate through the cells to place the next queen.
        for i in range(start_index, board_size * board_size):
            r = i // board_size
            c = i % board_size

            # Check if the current cell (r, c) is safe from attacks.
            is_safe = True
            for pr, pc in positions:
                # Check for row, column, and diagonal attacks (toroidal).
                if (r == pr or
                        c == pc or
                        (r - c) % board_size == (pr - pc) % board_size or
                        (r + c) % board_size == (pr + pc) % board_size):
                    is_safe = False
                    break
            
            if is_safe:
                # Place the queen and recurse.
                positions.append((r, c))
                backtrack(k + 1, i + 1, positions)
                # Backtrack: remove the queen to explore other possibilities.
                positions.pop()

    # Start the search from the first queen (k=0) at the first cell (start_index=0).
    backtrack(0, 0, [])

    print(f"Problem parameters:")
    print(f"Board size: {board_size}x{board_size} (toroidal)")
    print(f"Number of queens: {num_queens}")
    print(f"Final calculation: The number of ways is the result of the search.")
    print(f"Result: {count}")


solve_4_queens_toroidal()