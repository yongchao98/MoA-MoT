# Define the board state
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', 'X', 'X', 'X', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', 'O', 'O', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures stones
def check_capture(board, move):
    x, y = move
    # Simulate the move
    board[y][x] = 'X'
    # Check if the move captures any stones
    captured = False
    # Check the group at E9 (4, 8) and E8 (4, 7)
    if board[8][4] == 'O' and board[7][4] == 'O' and board[9][4] == 'X' and board[8][3] == 'X' and board[8][5] == 'X':
        captured = True
    return captured

# Check the move at F9 (5, 8)
move = (5, 8)
captured = check_capture(board, move)
print(captured)