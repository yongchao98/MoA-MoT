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

# Define the directions for adjacent positions
directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 11 and 0 <= y < 11

# Function to find liberties of a group
def find_liberties(x, y, visited):
    stack = [(x, y)]
    liberties = set()
    group = set()
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        group.add((cx, cy))
        for dx, dy in directions:
            nx, ny = cx + dx, cy + dy
            if is_within_board(nx, ny):
                if board[nx][ny] == '.':
                    liberties.add((nx, ny))
                elif board[nx][ny] == 'O' and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties, group

# Find all vulnerable white groups
vulnerable_groups = []
visited = set()
for i in range(11):
    for j in range(11):
        if board[i][j] == 'O' and (i, j) not in visited:
            liberties, group = find_liberties(i, j, visited)
            if len(liberties) == 1:
                vulnerable_groups.append((liberties.pop(), group))

# Determine the best move
best_move = None
max_capture = 0
for liberty, group in vulnerable_groups:
    if len(group) > max_capture:
        max_capture = len(group)
        best_move = liberty

# Convert the best move to board coordinates
if best_move:
    row, col = best_move
    move = f"{chr(col + ord('A'))}{11 - row}"
else:
    move = "No immediate capture"

print(move)