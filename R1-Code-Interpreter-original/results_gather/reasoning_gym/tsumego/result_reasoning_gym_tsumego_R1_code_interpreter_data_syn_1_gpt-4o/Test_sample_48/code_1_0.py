# Define the board
board = [
    ['X', '.', '.', '.', '.', 'X', 'X', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'O', 'X', '.', '.'],
    ['.', '.', '.', 'X', '.', 'X', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', 'X', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures any stones
def check_capture(board, move):
    x, y = move
    board[y][x] = 'X'  # Simulate the move
    captured = []

    # Check adjacent positions for potential captures
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(board[0]) and 0 <= ny < len(board):
            if board[ny][nx] == 'O':
                # Check if this group of 'O' stones is captured
                if is_captured(board, nx, ny):
                    captured.append((nx, ny))

    return captured

# Function to determine if a group of stones is captured
def is_captured(board, x, y):
    visited = set()
    stack = [(x, y)]
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < len(board[0]) and 0 <= ny < len(board):
                if board[ny][nx] == '.':
                    return False  # Found a liberty
                if board[ny][nx] == 'O' and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return True

# Check the potential moves
potential_moves = [(7, 9), (2, 5), (2, 1)]
captures = {move: check_capture(board, move) for move in potential_moves}

# Print the results
print(captures)