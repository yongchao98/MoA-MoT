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

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 11 and 0 <= y < 11

# Function to count liberties of a group
def count_liberties(x, y, color):
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

# Function to simulate a move and count captured stones
def simulate_move(x, y):
    if board[x][y] != '.':
        return 0
    # Place the stone
    board[x][y] = 'X'
    captured = 0
    # Check adjacent white groups
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nx, ny = x + dx, y + dy
        if is_within_board(nx, ny) and board[nx][ny] == 'O':
            if count_liberties(nx, ny, 'O') == 1:
                # Capture the group
                captured += 1
    # Remove the stone
    board[x][y] = '.'
    return captured

# Evaluate potential moves
potential_moves = [(6, 7), (5, 5), (7, 6), (4, 6)]
best_move = None
max_captured = 0
for move in potential_moves:
    captured = simulate_move(*move)
    if captured > max_captured:
        max_captured = captured
        best_move = move

# Convert the best move to board coordinates
if best_move:
    columns = "ABCDEFGHIJK"
    best_move_coord = f"{columns[best_move[1]]}{11 - best_move[0]}"
    print(best_move_coord)
else:
    print("No capturing move found")