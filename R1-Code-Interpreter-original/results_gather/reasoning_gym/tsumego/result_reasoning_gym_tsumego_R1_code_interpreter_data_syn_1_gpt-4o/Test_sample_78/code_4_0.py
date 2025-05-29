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
    # Check the group at D9 (3, 8) and E9 (4, 8)
    if board[8][3] == 'O' and board[8][4] == 'O' and board[7][3] == 'X' and board[9][3] == 'X' and board[8][2] == 'X' and board[8][5] == 'X':
        captured = True
    return captured

# Check the move at D8 (3, 7)
move = (3, 7)
captured = check_capture(board, move)
print(captured)