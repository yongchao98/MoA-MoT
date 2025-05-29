# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'O', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', 'X', 'O', 'X', '.', 'X', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures any stones
def check_capture(x, y, board):
    # Check if the move is valid
    if board[y][x] != '.':
        return 0
    
    # Simulate the move
    board[y][x] = 'X'
    
    # Check for captures
    captured = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 13 and 0 <= ny < 13 and board[ny][nx] == 'O':
            # Check if this group of 'O' is captured
            if is_captured(nx, ny, board):
                captured += count_group(nx, ny, board)
    
    # Undo the move
    board[y][x] = '.'
    return captured

# Function to check if a group is captured
def is_captured(x, y, board):
    visited = set()
    stack = [(x, y)]
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < 13 and 0 <= ny < 13:
                if board[ny][nx] == '.':
                    return False
                if board[ny][nx] == 'O' and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return True

# Function to count the number of stones in a group
def count_group(x, y, board):
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
            if 0 <= nx < 13 and 0 <= ny < 13 and board[ny][nx] == 'O' and (nx, ny) not in visited:
                stack.append((nx, ny))
    return count

# Check specific moves
moves_to_check = [(1, 12), (2, 11)]  # B1 and C2
best_move = None
max_captured = 0
for x, y in moves_to_check:
    captured = check_capture(x, y, board)
    if captured > max_captured:
        max_captured = captured
        best_move = (x, y)

# Convert the best move to Go coordinates
columns = "ABCDEFGHIJKLM"
if best_move:
    move = f"{columns[best_move[0]]}{13 - best_move[1]}"
else:
    move = "No move"

print(move)