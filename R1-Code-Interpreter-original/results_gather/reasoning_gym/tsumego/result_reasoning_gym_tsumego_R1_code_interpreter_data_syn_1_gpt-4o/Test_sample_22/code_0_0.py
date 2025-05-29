# Define the board
board = [
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', 'X', 'X', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'X', '.', '.', '.', '.', '.', '.', 'O'],
    ['X', 'X', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', 'X', 'O'],
    ['.', '.', '.', '.', '.', 'O', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 10 and 0 <= y < 10

# Function to get the liberties of a group
def get_liberties(x, y, color):
    visited = set()
    stack = [(x, y)]
    liberties = 0
    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if is_within_board(nx, ny):
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties

# Function to evaluate potential moves
def evaluate_moves():
    best_move = None
    max_captures = 0
    for x in range(10):
        for y in range(10):
            if board[x][y] == '.':
                # Simulate placing a black stone
                board[x][y] = 'X'
                captures = 0
                # Check each group of white stones
                for gx, gy in [(3, 9), (4, 8), (0, 5), (1, 5), (2, 5), (1, 4), (4, 3), (5, 1), (7, 1), (9, 4), (9, 2)]:
                    if board[gx][gy] == 'O' and get_liberties(gx, gy, 'O') == 0:
                        captures += 1
                # Revert the move
                board[x][y] = '.'
                # Update the best move
                if captures > max_captures:
                    max_captures = captures
                    best_move = (x, y)
    return best_move

# Find the best move
best_move = evaluate_moves()
# Convert to board coordinates
columns = 'ABCDEFGHIJ'
row, col = best_move
move = f"{columns[col]}{10 - row}"
print(move)