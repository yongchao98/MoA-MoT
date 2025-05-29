# Define the board
board = [
    ['.', '.', 'X', 'X', 'X', '.', 'O', '.', '.'],
    ['.', '.', 'X', 'O', 'O', 'X', '.', '.', 'X'],
    ['O', '.', 'X', 'O', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', '.', '.'],
    ['.', '.', 'X', '.', 'O', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures stones
def check_capture(board, move):
    x, y = move
    board[x][y] = 'X'  # Place the black stone
    captured = 0
    
    # Check surrounding positions for capture
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 9 and 0 <= ny < 9 and board[nx][ny] == 'O':
            # Check if the white stone is surrounded
            if all(0 <= nx + ddx < 9 and 0 <= ny + ddy < 9 and board[nx + ddx][ny + ddy] == 'X'
                   for ddx, ddy in directions if (ddx, ddy) != (-dx, -dy)):
                captured += 1
    
    board[x][y] = '.'  # Reset the board
    return captured

# Evaluate potential moves
moves = {'F7': (2, 5), 'E7': (2, 4), 'D6': (3, 3)}
captures = {move: check_capture(board, pos) for move, pos in moves.items()}

# Find the move with the maximum capture
best_move = max(captures, key=captures.get)
print(best_move)