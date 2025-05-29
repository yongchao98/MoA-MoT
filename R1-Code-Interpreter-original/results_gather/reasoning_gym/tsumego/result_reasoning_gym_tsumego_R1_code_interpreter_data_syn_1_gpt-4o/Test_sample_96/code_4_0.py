def simulate_move(board, move):
    # Simulate the move and check if it results in a capture
    x, y = move
    board[y][x] = 'X'  # Place the Black stone
    # Check for capture
    captured = False
    # Check surrounding White stones
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 10 and 0 <= ny < 10 and board[ny][nx] == 'O':
            # Check if this White stone is captured
            if is_captured(board, nx, ny):
                captured = True
    return captured

def is_captured(board, x, y):
    # Check if the White stone at (x, y) is captured
    visited = set()
    stack = [(x, y)]
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        # Check liberties
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < 10 and 0 <= ny < 10:
                if board[ny][nx] == '.':
                    return False  # Has a liberty
                elif board[ny][nx] == 'O':
                    stack.append((nx, ny))
    return True

# Initial board setup
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

# Test move
move_E8 = (4, 2)

capture_E8 = simulate_move([row[:] for row in board], move_E8)

print("Capture with E8:", capture_E8)