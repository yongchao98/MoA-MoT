# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.'],
    ['.', '.', '.', '.', 'X', 'O', 'O', '.', '.', '.'],
    ['.', '.', 'X', 'X', '.', 'X', 'O', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', 'X', 'O', '.', 'O']
]

# Function to check if a move captures stones
def is_capture_move(x, y, board):
    # Check if the move is within bounds and on an empty spot
    if 0 <= x < 10 and 0 <= y < 10 and board[y][x] == '.':
        # Simulate the move
        board[y][x] = 'X'
        # Check surrounding positions for capture
        captured = False
        directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 10 and 0 <= ny < 10 and board[ny][nx] == 'O':
                # Check if this white stone is surrounded
                if all(0 <= nx + ddx < 10 and 0 <= ny + ddy < 10 and (board[ny + ddy][nx + ddx] == 'X' or (nx + ddx == x and ny + ddy == y))
                       for ddx, ddy in directions if (ddx, ddy) != (-dx, -dy)):
                    captured = True
        # Revert the move
        board[y][x] = '.'
        return captured
    return False

# Check potential moves
potential_moves = [(x, y) for x in range(10) for y in range(10) if is_capture_move(x, y, board)]

# Print the potential moves
print(potential_moves)