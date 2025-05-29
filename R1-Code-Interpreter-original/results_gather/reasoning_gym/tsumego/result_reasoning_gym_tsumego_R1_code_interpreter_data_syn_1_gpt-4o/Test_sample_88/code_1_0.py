# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', 'O', 'X', 'X', '.'],
    ['.', 'X', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', 'X', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 13 and 0 <= y < 13

# Function to get adjacent positions
def get_adjacent(x, y):
    return [(x-1, y), (x+1, y), (x, y-1), (x, y+1)]

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
        for nx, ny in get_adjacent(cx, cy):
            if is_within_board(nx, ny):
                if board[nx][ny] == '.':
                    liberties += 1
                elif board[nx][ny] == color and (nx, ny) not in visited:
                    stack.append((nx, ny))
    return liberties

# Function to simulate a move and count captured stones
def simulate_move(x, y):
    if not is_within_board(x, y) or board[x][y] != '.':
        return 0
    # Place a black stone
    board[x][y] = 'X'
    captured = 0
    # Check adjacent white groups
    for nx, ny in get_adjacent(x, y):
        if is_within_board(nx, ny) and board[nx][ny] == 'O':
            if count_liberties(nx, ny, 'O') == 1:  # Check if placing a stone reduces liberties to zero
                # Capture the group
                captured += 1
    # Reset the board
    board[x][y] = '.'
    return captured

# Find the best move
best_move = None
max_captured = 0
for i in range(13):
    for j in range(13):
        captured = simulate_move(i, j)
        if captured > max_captured:
            max_captured = captured
            best_move = (i, j)

# Convert the best move to board coordinates
if best_move:
    row, col = best_move
    col_letter = chr(ord('A') + col)
    row_number = 13 - row
    print(f"<<<{col_letter}{row_number}>>>")
else:
    print("No capturing move found")