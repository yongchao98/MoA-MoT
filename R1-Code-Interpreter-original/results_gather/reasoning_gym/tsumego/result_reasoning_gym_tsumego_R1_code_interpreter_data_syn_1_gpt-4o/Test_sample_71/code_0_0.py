# Define the board
board = [
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'X'],
    ['.', 'X', '.', '.', '.', '.', 'X', '.', 'X', 'O', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O']
]

# Function to check if a move captures any stones
def check_capture(board, move):
    x, y = move
    board[x][y] = 'X'  # Place the black stone
    # Check if any white stones are captured
    captured = False
    # Check the group of white stones at the bottom
    white_group = [(10, 9), (10, 10), (11, 10), (11, 11)]
    liberties = 0
    for wx, wy in white_group:
        if board[wx][wy] == 'O':
            # Check liberties
            for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                nx, ny = wx + dx, wy + dy
                if 0 <= nx < 12 and 0 <= ny < 12 and board[nx][ny] == '.':
                    liberties += 1
    if liberties == 0:
        captured = True
    board[x][y] = '.'  # Reset the board
    return captured

# Check potential moves
moves = [(9, 9), (8, 10)]  # J3, K4
results = {move: check_capture(board, move) for move in moves}

# Print the results
print(results)