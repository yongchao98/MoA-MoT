# Define the board
board = [
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', 'X', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', 'O', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', 'X', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'X', 'O', 'O', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Directions: up, down, left, right
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Function to check if a move captures stones
def captures_stones(x, y, color):
    # Check if the move is within bounds
    if x < 0 or x >= len(board) or y < 0 or y >= len(board[0]):
        return 0
    # Check if the position is empty
    if board[x][y] != '.':
        return 0
    # Place the stone
    board[x][y] = color
    # Check for captures
    captured = 0
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(board) and 0 <= ny < len(board[0]) and board[nx][ny] == 'O':
            # Check if the group is surrounded
            if is_surrounded(nx, ny, 'O'):
                captured += capture_group(nx, ny)
    # Remove the stone
    board[x][y] = '.'
    return captured

# Function to check if a group is surrounded
def is_surrounded(x, y, color):
    visited = set()
    stack = [(x, y)]
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        for dx, dy in directions:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                if board[nx][ny] == '.':
                    return False
                if board[nx][ny] == color:
                    stack.append((nx, ny))
    return True

# Function to capture a group of stones
def capture_group(x, y):
    captured = 0
    stack = [(x, y)]
    while stack:
        cx, cy = stack.pop()
        if board[cx][cy] == 'O':
            board[cx][cy] = '.'
            captured += 1
            for dx, dy in directions:
                nx, ny = cx + dx, cy + dy
                if 0 <= nx < len(board) and 0 <= ny < len(board[0]):
                    stack.append((nx, ny))
    return captured

# Evaluate the specific move at C7
move = (5, 2)  # C7
captured = captures_stones(move[0], move[1], 'X')

# Convert the move to board coordinates
columns = 'ABCDEFGHIJKLM'
move_coord = f"{columns[move[1]]}{13 - move[0]}"
print(move_coord, captured)