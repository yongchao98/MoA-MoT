# Initial board setup
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', 'X', '.', '.', 'X', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', 'X', 'X', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', 'X', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', 'O', '.', '.', '.', '.', '.', 'O', '.', '.', 'O'],
    ['X', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.']
]

# Function to simulate the move and check capture
def simulate_move(board, move):
    x, y = move
    board[y][x] = 'X'  # Place the black stone
    # Check for capture around the move
    captured = []
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 12 and 0 <= ny < 12 and board[ny][nx] == 'O':
            # Check if this group of white stones is captured
            if is_captured(board, nx, ny):
                captured.append((nx, ny))
    return captured

def is_captured(board, x, y):
    # Simple check for capture: if no liberties, it's captured
    visited = set()
    stack = [(x, y)]
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < 12 and 0 <= ny < 12:
                if board[ny][nx] == '.':
                    return False  # Found a liberty
                if board[ny][nx] == 'O' and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return True

# Simulate the move at D4 (3, 4)
captured_stones = simulate_move(board, (3, 4))
print(captured_stones)