def find_optimal_moves():
    """
    Analyzes the given Connect 4 board to find all optimal moves for player 'O'
    to win as fast as possible.
    """
    board = [
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', '.', '.', '.', '.'],
        ['.', '.', '.', 'O', 'O', '.', '.'],
        ['O', '.', 'X', 'O', 'X', 'X', 'X'],
        ['X', 'O', 'O', 'X', 'X', 'O', 'X']
    ]
    player = 'O'
    rows = 6
    cols = 7
    col_names = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
    winning_moves = []

    def check_win(r, c):
        """Checks for a win condition around the given cell (r, c)."""
        # Check horizontal
        for j in range(cols - 3):
            if all(board[r][j+k] == player for k in range(4)):
                return True
        # Check vertical
        for i in range(rows - 3):
            if all(board[i+k][c] == player for k in range(4)):
                return True
        # Check positive diagonal (up-right)
        for i in range(rows - 3):
            for j in range(cols - 3):
                if all(board[i+k][j+k] == player for k in range(4)):
                    return True
        # Check negative diagonal (down-right)
        for i in range(3, rows):
            for j in range(cols - 3):
                if all(board[i-k][j+k] == player for k in range(4)):
                    return True
        return False

    # Iterate through each column to find possible moves
    for c in range(cols):
        # Find the lowest empty row in the current column
        r = -1
        for i in range(rows - 1, -1, -1):
            if board[i][c] == '.':
                r = i
                break
        
        # If column is not full, simulate the move
        if r != -1:
            # Temporarily place the piece
            board[r][c] = player
            
            # Check if this move results in a win
            if check_win(r, c):
                # Convert 0-indexed (r, c) to 1-indexed board notation (e.g., 'f3')
                move_name = f"{col_names[c]}{r + 1}"
                winning_moves.append(move_name)
            
            # Undo the move for the next iteration
            board[r][c] = '.'

    print(", ".join(winning_moves))

find_optimal_moves()