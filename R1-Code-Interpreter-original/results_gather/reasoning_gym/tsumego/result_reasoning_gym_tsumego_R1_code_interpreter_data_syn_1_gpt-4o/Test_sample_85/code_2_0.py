# Define the board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', 'X', '.', '.', '.', 'O', '.', 'X', '.'],
    ['.', '.', 'X', '.', '.', '.', 'X', 'X', '.', '.'],
    ['.', '.', '.', '.', 'O', 'X', 'O', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures any stones
def check_capture(board, move):
    x, y = move
    # Temporarily place the black stone
    board[x][y] = 'X'
    # Check the group of stones around the move
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    captured = False
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 10 and 0 <= ny < 10 and board[nx][ny] == 'O':
            # Check if this group of white stones is captured
            if is_captured(board, nx, ny):
                captured = True
    # Remove the temporary black stone
    board[x][y] = '.'
    return captured

# Function to check if a group of stones is captured
def is_captured(board, x, y):
    visited = set()
    stack = [(x, y)]
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        # Check for liberties
        directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        for dx, dy in directions:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < 10 and 0 <= ny < 10:
                if board[nx][ny] == '.':
                    return False
                elif board[nx][ny] == 'O' and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return True

# Check if playing at E5 captures any stones
move = (5, 4)  # F5
captured = check_capture(board, move)
print(captured)