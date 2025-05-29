# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', 'O', 'O', 'O', '.'],
    ['.', '.', '.', 'X', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', 'X', 'O', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', '.', '.', '.', '.', '.', 'X', '.', '.']
]

# Function to check if a move captures stones
def is_capture_move(board, x, y, player):
    opponent = 'O' if player == 'X' else 'X'
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    captured = False

    # Check if the move captures any opponent stones
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 10 and 0 <= ny < 10 and board[ny][nx] == opponent:
            # Check if the opponent group is captured
            if not has_liberty(board, nx, ny, opponent):
                captured = True
                break

    return captured

# Function to check if a group has any liberties
def has_liberty(board, x, y, player):
    visited = set()
    stack = [(x, y)]
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))

        for dx, dy in directions:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < 10 and 0 <= ny < 10:
                if board[ny][nx] == '.':
                    return True
                if board[ny][nx] == player and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return False

# Check if playing at E7 captures any stones
x, y = 4, 6  # E7
player = 'X'
capture = is_capture_move(board, x, y, player)
print(capture)