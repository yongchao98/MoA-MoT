# Define the board as a 2D list
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'O', 'O', 'X', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'X', '.', '.', 'O', '.', '.', '.', '.', '.']
]

# Function to check if a move captures any stones
def check_capture(x, y, color):
    # Directions: up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    captured = 0
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 13 and 0 <= ny < 13 and board[nx][ny] == 'O':
            # Check if the stone is surrounded
            surrounded = True
            for ddx, ddy in directions:
                nnx, nny = nx + ddx, ny + ddy
                if 0 <= nnx < 13 and 0 <= nny < 13 and board[nnx][nny] != 'X':
                    surrounded = False
                    break
            if surrounded:
                captured += 1
    return captured

# Check potential move at F1
capture_F1 = check_capture(12, 5, 'X')  # F1

# Print the result
print("Capture at F1:", capture_F1)