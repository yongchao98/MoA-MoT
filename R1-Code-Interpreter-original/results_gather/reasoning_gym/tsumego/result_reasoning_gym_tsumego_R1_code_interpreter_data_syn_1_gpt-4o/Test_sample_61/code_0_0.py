# Define the board state
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', 'X', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'O', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', 'O', '.', 'O', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures stones
def is_capture(move, board):
    # Simulate the move
    x, y = move
    board[y][x] = 'X'
    
    # Check if any white stones are captured
    # For simplicity, we will just check the immediate surroundings
    captured = False
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 12 and 0 <= ny < 12 and board[ny][nx] == 'O':
            # Check if this white stone is captured
            # For simplicity, assume it is captured if surrounded by black stones
            if all(0 <= nx + ddx < 12 and 0 <= ny + ddy < 12 and board[ny + ddy][nx + ddx] == 'X' for ddx, ddy in directions):
                captured = True
                break
    
    # Revert the move
    board[y][x] = '.'
    return captured

# Check potential moves
potential_moves = [(5, 5), (5, 7)]  # E6 and F8
for move in potential_moves:
    if is_capture(move, board):
        print(f"Move {move} captures stones.")