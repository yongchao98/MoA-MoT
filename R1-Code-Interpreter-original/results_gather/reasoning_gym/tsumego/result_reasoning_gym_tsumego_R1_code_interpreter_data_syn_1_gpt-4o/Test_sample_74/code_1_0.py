# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', 'X', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', 'X', 'O', 'O', 'X', '.', '.', 'X'],
    ['.', 'X', '.', '.', '.', 'X', 'O', '.', '.', '.', '.', '.'],
    ['.', 'X', 'X', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', 'X', '.']
]

# Function to check capture potential
def check_capture(x, y):
    # Check if the position is empty
    if board[x][y] != '.':
        return 0
    
    # Place a black stone
    board[x][y] = 'X'
    
    # Check for captures
    captured = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 12 and 0 <= ny < 12 and board[nx][ny] == 'O':
            # Check if this white stone is captured
            liberties = 0
            for ddx, ddy in directions:
                nnx, nny = nx + ddx, ny + ddy
                if 0 <= nnx < 12 and 0 <= nny < 12 and board[nnx][nny] == '.':
                    liberties += 1
            if liberties == 0:
                captured += 1
    
    # Remove the black stone
    board[x][y] = '.'
    
    return captured

# Evaluate potential moves
move_9C = check_capture(3, 2)
move_12J = check_capture(11, 9)

# Determine the best move
best_move = '<<<C9>>>' if move_9C > move_12J else '<<<J12>>>'
print(best_move)