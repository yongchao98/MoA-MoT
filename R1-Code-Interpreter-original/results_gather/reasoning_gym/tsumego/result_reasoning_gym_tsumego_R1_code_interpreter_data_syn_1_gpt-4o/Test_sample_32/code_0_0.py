# Define the board as a list of lists
board = [
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', 'X', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O', 'X'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', 'X', 'O', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check liberties of a group
def check_liberties(x, y, color):
    visited = set()
    stack = [(x, y)]
    liberties = 0
    group = []

    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        group.append((cx, cy))

        # Check all 4 directions
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if 0 <= nx < 12 and 0 <= ny < 12:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return liberties, group

# Find the best move for Black
best_move = None
max_captures = 0

for i in range(12):
    for j in range(12):
        if board[i][j] == 'O':
            liberties, group = check_liberties(i, j, 'O')
            if liberties == 1:  # If the group has only one liberty, it can be captured
                for x, y in group:
                    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                        nx, ny = x + dx, y + dy
                        if 0 <= nx < 12 and 0 <= ny < 12 and board[nx][ny] == '.':
                            if len(group) > max_captures:
                                max_captures = len(group)
                                best_move = (nx, ny)

# Convert the best move to board coordinates
if best_move:
    row, col = best_move
    move = f"{chr(col + ord('A'))}{12 - row}"
else:
    move = "No immediate capture possible"

print(move)