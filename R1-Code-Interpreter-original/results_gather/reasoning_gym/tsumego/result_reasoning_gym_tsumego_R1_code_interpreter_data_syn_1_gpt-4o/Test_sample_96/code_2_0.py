# Define the board
board = [
    ['.', '.', 'X', '.', 'X', '.', '.', '.', '.', '.'],
    ['X', 'X', 'O', 'O', 'O', 'X', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', 'O', 'X', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['O', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.']
]

# Function to check if a stone has liberties
def has_liberties(board, x, y):
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 10 and 0 <= ny < 10 and board[ny][nx] == '.':
            return True
    return False

# Function to simulate the move and check captures
def simulate_move(board, move):
    x, y = move
    board[y][x] = 'X'  # Place the black stone
    captures = []
    # Check the surrounding positions
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 10 and 0 <= ny < 10 and board[ny][nx] == 'O':
            # Check if the white stone is captured
            if not has_liberties(board, nx, ny):
                captures.append((nx, ny))
    return captures

# Simulate the move at C8
move = (2, 2)  # C8 in 0-indexed (x, y) format
captures = simulate_move(board, move)
print("Captures:", captures)