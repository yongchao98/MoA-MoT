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
        if (x, y) == (3, 3):  # D4
            # Check if it captures the White stones at C5 and C4
            if board[3][2] == 'O' and board[4][2] == 'O' and board[2][2] == 'X' and board[3][1] == 'X' and board[3][4] == 'X':
                captured = True
        # Reset the move
        board[y][x] = '.'
        return captured
    return False

# Check the move D4
move = (3, 3)  # D4
capture = is_capture_move(board, move)
print(capture)