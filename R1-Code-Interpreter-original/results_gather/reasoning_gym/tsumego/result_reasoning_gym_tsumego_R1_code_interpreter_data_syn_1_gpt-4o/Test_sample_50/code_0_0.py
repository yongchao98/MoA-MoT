# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', 'O', '.', '.', '.', '.'],
    ['X', 'O', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', 'O', 'X', '.', '.', 'O', '.', '.', '.'],
    ['X', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to find groups of stones
def find_groups(board, stone):
    visited = set()
    groups = []

    def dfs(x, y, group):
        if (x, y) in visited or board[y][x] != stone:
            return
        visited.add((x, y))
        group.append((x, y))
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 10 and 0 <= ny < 10:
                dfs(nx, ny, group)

    for y in range(10):
        for x in range(10):
            if (x, y) not in visited and board[y][x] == stone:
                group = []
                dfs(x, y, group)
                groups.append(group)

    return groups

# Function to count liberties of a group
def count_liberties(group):
    liberties = set()
    for x, y in group:
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 10 and 0 <= ny < 10 and board[ny][nx] == '.':
                liberties.add((nx, ny))
    return liberties

# Find all white groups and their liberties
white_groups = find_groups(board, 'O')
best_move = None
max_captured = 0

for group in white_groups:
    liberties = count_liberties(group)
    if len(liberties) == 1:  # Group is capturable
        liberty = liberties.pop()
        if len(group) > max_captured:
            max_captured = len(group)
            best_move = liberty

# Convert best move to board coordinates
if best_move:
    x, y = best_move
    move = f"{chr(x + ord('A'))}{10 - y}"
else:
    move = "No immediate capture"

print(move)