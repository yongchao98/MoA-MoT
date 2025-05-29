# Define the board
board = [
    ['.', '.', '.', '.', 'X', '.', '.', 'X', 'X'],
    ['.', '.', '.', '.', '.', 'X', 'X', 'O', 'O'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', '.'],
    ['.', 'O', 'O', '.', '.', '.', 'X', 'X', 'O'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.']
]

# Function to check if a move at (x, y) captures any stones
def check_capture(x, y):
    # Directions: up, down, left, right
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    captured_stones = 0
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 9 and 0 <= ny < 9 and board[nx][ny] == 'O':
            # Check if this group of 'O' is captured
            if is_captured(nx, ny):
                captured_stones += count_group(nx, ny)
    return captured_stones

# Function to check if a group is captured
def is_captured(x, y):
    visited = set()
    return not has_liberty(x, y, visited)

# Function to check if a group has liberty
def has_liberty(x, y, visited):
    if (x, y) in visited:
        return False
    visited.add((x, y))
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 9 and 0 <= ny < 9:
            if board[nx][ny] == '.':
                return True
            if board[nx][ny] == 'O' and has_liberty(nx, ny, visited):
                return True
    return False

# Function to count the number of stones in a group
def count_group(x, y):
    visited = set()
    return count_stones(x, y, visited)

def count_stones(x, y, visited):
    if (x, y) in visited:
        return 0
    visited.add((x, y))
    count = 1
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 9 and 0 <= ny < 9 and board[nx][ny] == 'O':
            count += count_stones(nx, ny, visited)
    return count

# Check the move at G8
capture_count = check_capture(7, 6)
print(capture_count)