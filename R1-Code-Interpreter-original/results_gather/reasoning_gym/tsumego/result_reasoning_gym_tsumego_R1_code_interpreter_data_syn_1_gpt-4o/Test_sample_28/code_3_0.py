# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', 'X', '.', '.', 'X', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', 'O', 'X', 'X', '.', '.', '.'],
    ['.', '.', 'X', 'X', 'O', 'O', '.', '.', '.', 'X', '.'],
    ['O', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X']
]

# Function to check if a move captures any stones
def check_capture(x, y, color):
    # Directions: up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    captured = 0
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 11 and 0 <= ny < 11 and board[nx][ny] == 'O':
            # Check if this group of 'O' is captured
            if is_captured(nx, ny):
                captured += 1
    return captured

# Function to determine if a group is captured
def is_captured(x, y):
    # Use BFS to check if the group is surrounded
    queue = [(x, y)]
    visited = set(queue)
    while queue:
        cx, cy = queue.pop(0)
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < 11 and 0 <= ny < 11:
                if board[nx][ny] == '.':
                    return False
                if board[nx][ny] == 'O' and (nx, ny) not in visited:
                    visited.add((nx, ny))
                    queue.append((nx, ny))
    return True

# Check the move at G7
move_x, move_y = 4, 6  # G7
captured_stones = check_capture(move_x, move_y, 'X')
print(captured_stones)