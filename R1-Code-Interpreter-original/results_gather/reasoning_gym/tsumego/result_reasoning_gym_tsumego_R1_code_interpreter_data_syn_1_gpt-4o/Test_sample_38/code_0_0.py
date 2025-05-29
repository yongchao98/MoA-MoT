# Define the board
board = [
    ['.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['X', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', 'O', '.', '.'],
    ['.', '.', 'X', 'O', 'X', '.', '.', 'X', '.'],
    ['.', 'X', 'O', 'O', '.', 'O', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'X', '.', '.', 'O', '.'],
    ['.', 'X', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', 'O', '.', '.', '.']
]

# Function to find groups and their liberties
def find_groups_and_liberties(board):
    def get_neighbors(r, c):
        neighbors = []
        if r > 0: neighbors.append((r-1, c))
        if r < 8: neighbors.append((r+1, c))
        if c > 0: neighbors.append((r, c-1))
        if c < 8: neighbors.append((r, c+1))
        return neighbors

    visited = set()
    groups = []

    for r in range(9):
        for c in range(9):
            if board[r][c] == 'O' and (r, c) not in visited:
                # Perform a BFS to find the group
                queue = [(r, c)]
                group = []
                liberties = set()
                while queue:
                    cr, cc = queue.pop(0)
                    if (cr, cc) in visited:
                        continue
                    visited.add((cr, cc))
                    group.append((cr, cc))
                    for nr, nc in get_neighbors(cr, cc):
                        if board[nr][nc] == 'O' and (nr, nc) not in visited:
                            queue.append((nr, nc))
                        elif board[nr][nc] == '.':
                            liberties.add((nr, nc))
                groups.append((group, liberties))
    return groups

# Find groups and their liberties
groups = find_groups_and_liberties(board)

# Evaluate potential moves
best_move = None
max_captures = 0

for group, liberties in groups:
    if len(liberties) == 1:
        # If there's only one liberty, we can capture the group
        liberty = liberties.pop()
        if len(group) > max_captures:
            max_captures = len(group)
            best_move = liberty

# Convert the best move to the board coordinates
if best_move:
    row, col = best_move
    move = chr(col + ord('A')) + str(9 - row)
else:
    move = None

print(move)