# Define the board
board = [
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', 'X', '.', '.', 'X', '.'],
    ['.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.'],
    ['O', '.', 'X', 'O', 'X', '.', '.', '.', '.'],
    ['O', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', 'O', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', 'X', '.', '.', '.']
]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 9 and 0 <= y < 9

# Function to get the liberties of a group
def get_liberties(x, y, visited):
    stack = [(x, y)]
    liberties = set()
    group = set()
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        group.add((cx, cy))
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if is_within_board(nx, ny):
                if board[nx][ny] == '.':
                    liberties.add((nx, ny))
                elif board[nx][ny] == 'O' and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties, group

# Find all groups of white stones and their liberties
white_groups = []
visited = set()
for i in range(9):
    for j in range(9):
        if board[i][j] == 'O' and (i, j) not in visited:
            liberties, group = get_liberties(i, j, visited)
            white_groups.append((liberties, group))

# Determine the best move for Black
best_move = None
max_captured = 0
for liberties, group in white_groups:
    if len(liberties) == 1:  # Only one liberty means the group can be captured
        liberty = liberties.pop()
        if len(group) > max_captured:
            max_captured = len(group)
            best_move = liberty

# Convert the best move to board coordinates
if best_move:
    row, col = best_move
    move = f"{chr(col + ord('A'))}{9 - row}"
else:
    move = "No capturing move found"

print(move)