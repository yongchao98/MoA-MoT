# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', 'X', '.', '.', 'X', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', 'O', 'X', 'X', '.', '.', '.'],
    ['.', '.', 'X', 'X', 'O', 'O', '.', '.', '.', 'X', '.'],
    ['O', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X']
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

        # Check adjacent positions
        for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
            if 0 <= nx < 11 and 0 <= ny < 11:
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))

    return liberties, group

# Find the best move
best_move = None
max_captured = 0

for i in range(11):
    for j in range(11):
        if board[i][j] == 'O':
            liberties, group = check_liberties(i, j, 'O')
            if liberties == 1:
                # Find the liberty position
                for x, y in group:
                    for nx, ny in [(x-1, y), (x+1, y), (x, y-1), (x, y+1)]:
                        if 0 <= nx < 11 and 0 <= ny < 11 and board[nx][ny] == '.':
                            if len(group) > max_captured:
                                max_captured = len(group)
                                best_move = (nx, ny)

# Convert best move to board coordinates
if best_move:
    row = 11 - best_move[0]
    col = chr(ord('A') + best_move[1])
    print(f'<<<{col}{row}>>>')
else:
    print("No immediate capture possible")