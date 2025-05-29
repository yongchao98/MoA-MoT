# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', 'X', 'O', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', 'X', 'O', 'X'],
    ['.', '.', 'O', 'X', 'O', 'O', 'O', 'X', '.'],
    ['.', '.', 'X', '.', 'X', 'O', 'X', '.', '.'],
    ['.', 'O', 'X', 'X', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures any stones
def check_capture(board, move):
    x, y = move
    # Simulate the move
    board[y][x] = 'X'
    # Check for captures
    captures = 0
    # Check adjacent positions for potential captures
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 9 and 0 <= ny < 9 and board[ny][nx] == 'O':
            # Check if the group of stones is captured
            if is_captured(board, nx, ny):
                captures += 1
    # Revert the move
    board[y][x] = '.'
    return captures

# Function to determine if a group of stones is captured
def is_captured(board, x, y):
    visited = set()
    stack = [(x, y)]
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        # Check for liberties
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < 9 and 0 <= ny < 9:
                if board[ny][nx] == '.':
                    return False
                elif board[ny][nx] == 'O':
                    stack.append((nx, ny))
    return True

# Evaluate potential moves
potential_moves = {'D7': (3, 2), 'E7': (4, 2), 'F8': (5, 1)}
captures = {move: check_capture(board, coord) for move, coord in potential_moves.items()}

# Find the move with the maximum captures
best_move = max(captures, key=captures.get)
print(best_move)