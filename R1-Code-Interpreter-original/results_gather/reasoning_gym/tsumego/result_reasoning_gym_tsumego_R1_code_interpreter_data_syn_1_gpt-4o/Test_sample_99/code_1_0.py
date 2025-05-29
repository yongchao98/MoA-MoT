# Define the board
board = [
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', 'X', 'O', 'X', 'O', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'O', '.', 'X', '.', 'O', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'X', '.', '.', '.', 'X', '.', 'X'],
    ['X', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures any stones
def check_capture(board, move):
    x, y = move
    board[x][y] = 'X'  # Place the black stone
    # Check for captures around the move
    captures = 0
    # Directions: up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(board) and 0 <= ny < len(board[0]) and board[nx][ny] == 'O':
            # Check if this group of white stones is captured
            if is_captured(board, nx, ny):
                captures += count_group(board, nx, ny)
    return captures

# Function to check if a group of stones is captured
def is_captured(board, x, y):
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
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    return False  # Has a liberty
                elif board[nx][ny] == 'O':
                    stack.append((nx, ny))
    return True

# Function to count the number of stones in a group
def count_group(board, x, y):
    visited = set()
    stack = [(x, y)]
    count = 0
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        count += 1
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]) and board[nx][ny] == 'O':
                stack.append((nx, ny))
    return count

# Check the move at C10
move = (1, 1)  # C10
captures = check_capture(board, move)
print(captures)