# Define the board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', '.', '.', '.', '.', '.', '.'],
    ['X', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', 'O', 'X', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', 'X', '.', '.', '.', '.', 'X', '.', '.'],
    ['O', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.']
]

# Function to check if a group of stones is captured
def is_captured(board, x, y, color):
    visited = set()
    stack = [(x, y)]
    captured = True
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        # Check if the stone has a liberty
        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < 9 and 0 <= ny < 9:
                if board[ny][nx] == '.':
                    captured = False
                elif board[ny][nx] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return captured

# Function to simulate the move and check captures
def simulate_move(board, move):
    x, y = move
    board[y][x] = 'X'  # Place the black stone
    captures = 0
    # Check all white stones to see if they are captured
    for i in range(9):
        for j in range(9):
            if board[j][i] == 'O' and is_captured(board, i, j, 'O'):
                captures += 1
                board[j][i] = '.'  # Remove captured stone
    return captures

# Simulate the move at D5 (3, 4 in 0-indexed)
captures = simulate_move(board, (3, 4))
print(captures)