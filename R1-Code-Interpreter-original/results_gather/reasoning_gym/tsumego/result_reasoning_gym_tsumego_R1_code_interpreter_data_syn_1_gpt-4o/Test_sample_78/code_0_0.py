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
    # Check the group at E8 (4, 7)
    if board[7][4] == 'O' and board[6][4] == 'X' and board[8][4] == 'X' and board[7][3] == 'X' and board[7][5] == 'X':
        captured = True
    return captured

# Check the move at E7 (4, 7)
move = (4, 7)
captured = check_capture(board, move)
print(captured)