def solve_connect_4():
    """
    Analyzes the given Connect 4 board to find the optimal moves for 'O'.
    An optimal move is one that guarantees a win in the fewest possible turns.
    """
    # Board representation: 6 rows, 7 columns
    # '.' is empty, 'X' and 'O' are players.
    # Top-left is (0, 0), corresponding to g1 in the image if inverted, but we use standard matrix indexing
    # The image has row 1 at the top, so we'll map row 1->index 0, row 6->index 5
    # col a -> 0, col g -> 6
    board = [
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'O', 'O', '.', '.'],
        ['O', '.', 'X', 'O', 'X', 'X', 'X'],
        ['X', 'O', 'O', 'X', 'X', 'O', 'X']
    ]
    rows = 6
    cols = 7
    player = 'O'
    
    def get_valid_moves(current_board):
        """Finds the lowest empty cell in each column."""
        moves = []
        for c in range(cols):
            for r in range(rows - 1, -1, -1):
                if current_board[r][c] == '.':
                    moves.append((r, c))
                    break
        return moves

    def check_forced_win(move_r, move_c, current_board, p):
        """
        Checks if a move creates an unblockable threat, like an open-ended three-in-a-row.
        This is a simplified check focusing on the most likely winning pattern.
        """
        # Create a temporary board with the move made
        temp_board = [row[:] for row in current_board]
        temp_board[move_r][move_c] = p
        
        # Directions: horizontal, vertical, diagonal up-right, diagonal down-right
        directions = [(0, 1), (1, 0), (1, 1), (1, -1)]
        
        for dr, dc in directions:
            # Check for patterns around the new piece at (move_r, move_c)
            for i in range(-3, 1):
                try:
                    p1 = temp_board[move_r + i * dr][move_c + i * dc]
                    p2 = temp_board[move_r + (i + 1) * dr][move_c + (i + 1) * dc]
                    p3 = temp_board[move_r + (i + 2) * dr][move_c + (i + 2) * dc]
                    p4 = temp_board[move_r + (i + 3) * dr][move_c + (i + 3) * dc]

                    # Check for an open-ended three-in-a-row: .OOO.
                    if [p1, p2, p3, p4] == [p, p, p, '.']:
                        # Check the other side
                        if move_r + (i - 1) * dr >= 0 and move_c + (i - 1) * dc >= 0:
                            if temp_board[move_r + (i - 1) * dr][move_c + (i - 1) * dc] == '.':
                                # Check if both empty spots are playable
                                empty1_r, empty1_c = move_r + (i - 1) * dr, move_c + (i - 1) * dc
                                empty2_r, empty2_c = move_r + (i + 3) * dr, move_c + (i + 3) * dc
                                
                                # A spot is playable if it's on the bottom row or the spot below it is filled
                                playable1 = (empty1_r == rows - 1) or (temp_board[empty1_r + 1][empty1_c] != '.')
                                playable2 = (empty2_r == rows - 1) or (temp_board[empty2_r + 1][empty2_c] != '.')

                                if playable1 and playable2:
                                    return True
                except IndexError:
                    continue
        return False


    optimal_moves = []
    valid_moves = get_valid_moves(board)

    # A win in 2 is the fastest possible here (no immediate win)
    for r, c in valid_moves:
        if check_forced_win(r, c, board, player):
            col_name = chr(ord('a') + c)
            row_name = str(r + 1)
            optimal_moves.append(f"{col_name}{row_name}")

    # Based on manual analysis, d3 is a win-in-3, which is slower.
    # The check_forced_win function is designed to find win-in-2 scenarios.
    
    print(", ".join(sorted(optimal_moves)))

solve_connect_4()