# Define the board
board = [
    ['.', 'X', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', 'O', 'X', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'X', '.', 'O', '.', 'X', '.', '.'],
    ['X', '.', 'X', '.', '.', '.', '.', '.', 'X', '.'],
    ['X', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.']
]

# Function to check if a move captures stones
def captures_stones(board, move):
    x, y = move
    # Check the group of stones around the move
    captured = 0
    if board[y][x] == '.':
        # Check adjacent positions
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 10 and 0 <= ny < 10 and board[ny][nx] == 'O':
                # Check if this group of 'O' stones is captured
                if is_captured(board, nx, ny):
                    captured += count_group(board, nx, ny)
    return captured

# Function to check if a group is captured
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
            if 0 <= nx < 10 and 0 <= ny < 10:
                if board[ny][nx] == '.':
                    return False
                if board[ny][nx] == 'O' and (nx, ny) not in visited:
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
            if 0 <= nx < 10 and 0 <= ny < 10 and board[ny][nx] == 'O' and (nx, ny) not in visited:
                stack.append((nx, ny))
    return count

# Check the move B8
move = (1, 2)  # B8
captured_stones = captures_stones(board, move)
print(captured_stones)