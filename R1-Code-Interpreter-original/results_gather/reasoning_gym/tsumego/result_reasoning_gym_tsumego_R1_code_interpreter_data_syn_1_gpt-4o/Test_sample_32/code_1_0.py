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

# Find a strategic move for Black
best_move = None
min_liberties = float('inf')

for i in range(12):
    for j in range(12):
        if board[i][j] == 'O':
            liberties, group = check_liberties(i, j, 'O')
            if liberties < min_liberties:
                for x, y in group:
                    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                        nx, ny = x + dx, y + dy
                        if 0 <= nx < 12 and 0 <= ny < 12 and board[nx][ny] == '.':
                            min_liberties = liberties
                            best_move = (nx, ny)

# Convert the best move to board coordinates
if best_move:
    row, col = best_move
    move = f"{chr(col + ord('A'))}{12 - row}"
else:
    move = "No strategic move found"

print(move)