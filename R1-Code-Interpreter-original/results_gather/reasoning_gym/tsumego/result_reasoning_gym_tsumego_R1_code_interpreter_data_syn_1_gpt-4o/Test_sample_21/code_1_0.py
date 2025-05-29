# Define the board
board = [
    ['.', '.', 'O', '.', '.', 'X', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O', 'O'],
    ['.', '.', 'X', '.', '.', '.', '.', 'O', '.', '.', 'X', 'O', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures stones
def check_capture(x, y):
    # Check if the move at (x, y) captures any stones
    # For simplicity, we assume the move is valid and on an empty spot
    # Check the surrounding stones and their liberties
    captured = False
    if board[y][x] == '.':
        # Place the stone
        board[y][x] = 'X'
        # Check the surrounding positions
        directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]
        for dx, dy in directions:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 13 and 0 <= ny < 13 and board[ny][nx] == 'O':
                # Check if the white group is captured
                if is_captured(nx, ny):
                    captured = True
        # Remove the stone after checking
        board[y][x] = '.'
    return captured

# Function to check if a group is captured
def is_captured(x, y):
    # Use a simple flood fill to check if the group has no liberties
    visited = set()
    stack = [(x, y)]
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        # Check liberties
        directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]
        for dx, dy in directions:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < 13 and 0 <= ny < 13:
                if board[ny][nx] == '.':
                    return False
                if board[ny][nx] == 'O' and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return True

# Check the move at H10
move_x, move_y = 7, 9  # H10
captured = check_capture(move_x, move_y)
print(captured)