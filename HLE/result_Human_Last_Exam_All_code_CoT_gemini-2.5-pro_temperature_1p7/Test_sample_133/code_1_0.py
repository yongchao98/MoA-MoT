def find_optimal_moves():
    """
    Analyzes the given Connect 4 board to find the fastest winning moves for 'O'.
    It assumes pieces can be placed in any empty spot (no gravity).
    """

    # The board is represented as a list of lists.
    # board[row_index][col_index] corresponds to the image.
    # Row 1 is index 0, Row 6 is index 5.
    # Column 'a' is index 0, Column 'g' is index 6.
    board = [
    #   a    b    c    d    e    f    g
        ['.', '.', '.', '.', '.', '.', '.'], # Row 1
        ['.', '.', '.', '.', '.', '.', '.'], # Row 2
        ['.', '.', '.', '.', '.', '.', '.'], # Row 3
        ['.', '.', '.', 'O', 'O', '.', '.'], # Row 4
        ['O', '.', 'X', 'O', 'X', 'X', 'X'], # Row 5
        ['X', 'O', 'O', 'X', 'X', 'O', 'X']  # Row 6
    ]

    player = 'O'
    rows = 6
    cols = 7
    winning_moves = []

    def check_win_from_move(b, p, r, c):
        """
        Checks if a move for player 'p' at (r, c) on board 'b' results in a win.
        """
        # Directions: horizontal, vertical, diagonal, anti-diagonal
        directions = [(0, 1), (1, 0), (1, 1), (1, -1)]
        for dr, dc in directions:
            count = 1
            # Count in the positive direction
            for i in range(1, 4):
                nr, nc = r + i * dr, c + i * dc
                if 0 <= nr < rows and 0 <= nc < cols and b[nr][nc] == p:
                    count += 1
                else:
                    break
            # Count in the negative direction
            for i in range(1, 4):
                nr, nc = r - i * dr, c - i * dc
                if 0 <= nr < rows and 0 <= nc < cols and b[nr][nc] == p:
                    count += 1
                else:
                    break
            if count >= 4:
                return True
        return False

    # Find all empty cells on the board
    empty_cells = []
    for r in range(rows):
        for c in range(cols):
            if board[r][c] == '.':
                empty_cells.append((r, c))

    # Test each empty cell for a potential win in one move
    for r, c in empty_cells:
        temp_board = [row[:] for row in board]
        temp_board[r][c] = player
        
        if check_win_from_move(temp_board, player, r, c):
            # Convert 0-indexed (r, c) to "a1" format
            col_char = chr(ord('a') + c)
            row_char = str(r + 1)
            winning_moves.append(f"{col_char}{row_char}")
            
    # The fastest possible win is a win in one move.
    if winning_moves:
        # Sort for consistent output and join with commas
        print(", ".join(sorted(winning_moves)))
    else:
        # This case is not reached for the given board, as a win-in-1 exists.
        # A more complex analysis for forced wins would be required otherwise.
        print("No winning move in one turn found.")

find_optimal_moves()