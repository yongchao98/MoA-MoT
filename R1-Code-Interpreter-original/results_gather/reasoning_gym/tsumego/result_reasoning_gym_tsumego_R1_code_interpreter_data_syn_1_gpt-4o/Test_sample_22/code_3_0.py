# Define the board
board = [
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', 'X', 'X', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'X', '.', '.', '.', '.', '.', '.', 'O'],
    ['X', 'X', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', 'X', 'O'],
    ['.', '.', '.', '.', '.', 'O', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move at (x, y) captures the group at (gx, gy)
def check_capture(x, y, gx, gy):
    # Place the stone
    board[x][y] = 'X'
    
    # Check if the group at (gx, gy) is captured
    if is_captured(gx, gy, 'O'):
        return count_group(gx, gy, 'O')
    return 0

# Function to determine if a group is captured
def is_captured(x, y, color):
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
                if board[nx][ny] == '.':
                    return False
                if board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return True

# Function to count the number of stones in a group
def count_group(x, y, color):
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
            if 0 <= nx < 10 and 0 <= ny < 10 and board[nx][ny] == color and (nx, ny) not in visited:
                stack.append((nx, ny))
    return count

# Check the move at C5 (row 5, column 3) for capturing the group at A6 (row 6, column 1)
move_x, move_y = 4, 2
group_x, group_y = 4, 0
captured_stones = check_capture(move_x, move_y, group_x, group_y)
print(captured_stones)