# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['O', '.', 'X', 'O', 'X', '.', '.', '.', 'X'],
    ['O', 'X', 'O', 'O', '.', '.', '.', '.', '.'],
    ['.', 'X', 'X', 'O', 'X', 'X', '.', '.', '.'],
    ['X', 'X', '.', 'X', '.', 'O', '.', '.', '.']
]

# Function to check if a move captures stones
def check_capture(x, y, color):
    # Directions: up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    opponent = 'O' if color == 'X' else 'X'
    captured = False
    
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 9 and 0 <= ny < 9 and board[ny][nx] == opponent:
            # Check if the opponent's stone is in atari
            liberties = 0
            for ddx, ddy in directions:
                nnx, nny = nx + ddx, ny + ddy
                if 0 <= nnx < 9 and 0 <= nny < 9 and board[nny][nnx] == '.':
                    liberties += 1
            if liberties == 0:
                captured = True
    return captured

# Check the move at E2
move = 'E2'
x, y = 4, 1
result = check_capture(x, y, 'X')

# Print the result
print({move: result})