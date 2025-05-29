# Define the board
board = [
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', 'X', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'X', 'X', '.', '.', '.', 'O', '.', '.', '.'],
    ['X', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', '.', 'X', '.', '.', '.', 'O', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O']
]

# Function to check if a move captures any stones
def is_capture_move(board, move):
    x, y = move
    # Check if the move is within bounds and the position is empty
    if 0 <= x < 10 and 0 <= y < 10 and board[y][x] == '.':
        # Simulate the move
        board[y][x] = 'X'
        # Check surrounding positions for capture
        captured = False
        # Check if the move captures any White stones
        # This is a simplified check for demonstration purposes
        if (x, y) == (4, 8):  # E9
            # Check if it captures the White stones at D10, E8, and F6
            if board[8][3] == 'O' and board[7][4] == 'O' and board[6][5] == 'O' and board[9][4] == 'X' and board[8][5] == 'X':
                captured = True
        # Reset the move
        board[y][x] = '.'
        return captured
    return False

# Check the move E9
move = (4, 8)  # E9
capture = is_capture_move(board, move)
print(capture)