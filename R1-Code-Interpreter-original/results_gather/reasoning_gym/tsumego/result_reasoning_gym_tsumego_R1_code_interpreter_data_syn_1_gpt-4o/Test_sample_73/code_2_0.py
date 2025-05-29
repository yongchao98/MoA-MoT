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

# Evaluate potential moves
moves = [(5, 2), (10, 10), (3, 5), (3, 6)]
best_move = None
max_captured = 0
for x, y in moves:
    captured = captures_stones(x, y, 'X')
    if captured > max_captured:
        max_captured = captured
        best_move = (x, y)

# Convert the best move to board coordinates
if best_move is not None:
    columns = 'ABCDEFGHIJKLM'
    best_move_coord = f"{columns[best_move[1]]}{13 - best_move[0]}"
    print(best_move_coord)
else:
    print("No capturing move found")