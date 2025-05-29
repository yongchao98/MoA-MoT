# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'X', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['O', 'O', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'X', '.', '.', 'O', '.', '.', '.', 'X', '.'],
    ['.', 'X', '.', '.', '.', '.', 'O', '.', '.', '.', 'X'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to find groups and their liberties
def find_groups_and_liberties(board):
    def get_neighbors(x, y):
        neighbors = []
        if x > 0: neighbors.append((x-1, y))
        if x < len(board) - 1: neighbors.append((x+1, y))
        if y > 0: neighbors.append((x, y-1))
        if y < len(board[0]) - 1: neighbors.append((x, y+1))
        return neighbors

    visited = set()
    groups = []

    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] == 'O' and (i, j) not in visited:
                stack = [(i, j)]
                group = []
                liberties = set()
                while stack:
                    x, y = stack.pop()
                    if (x, y) in visited:
                        continue
                    visited.add((x, y))
                    group.append((x, y))
                    for nx, ny in get_neighbors(x, y):
                        if board[nx][ny] == 'O' and (nx, ny) not in visited:
                            stack.append((nx, ny))
                        elif board[nx][ny] == '.':
                            liberties.add((nx, ny))
                groups.append((group, liberties))
    return groups

# Find groups and their liberties
groups = find_groups_and_liberties(board)

# Determine the best move
best_move = None
max_captures = 0

for group, liberties in groups:
    if len(liberties) == 1:
        # Group is in atari, can be captured
        liberty = liberties.pop()
        if len(group) > max_captures:
            max_captures = len(group)
            best_move = liberty

# Convert best move to board coordinates
if best_move:
    row, col = best_move
    move = f"<<<{chr(col + ord('A'))}{11 - row}>>>"
else:
    move = "No immediate capture possible"

print(move)