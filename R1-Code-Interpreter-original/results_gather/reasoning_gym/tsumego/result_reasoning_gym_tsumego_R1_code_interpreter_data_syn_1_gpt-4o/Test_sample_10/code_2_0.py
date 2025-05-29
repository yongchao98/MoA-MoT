# Define the board
board = [
    ['X', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', 'O', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O'],
    ['O', '.', 'X', 'O', 'X', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'O', '.', 'O', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', 'O', '.', '.', 'X', '.', '.']
]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 9 and 0 <= y < 9

# Function to check if a group of stones is captured
def is_captured(x, y, color):
    visited = set()
    stack = [(x, y)]
    liberties = 0

    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))

        # Check the four adjacent positions
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = cx + dx, cy + dy
            if is_within_board(nx, ny):
                if board[ny][nx] == '.':
                    liberties += 1
                elif board[ny][nx] == color:
                    stack.append((nx, ny))

    return liberties == 0

# Check potential moves
potential_moves = [(4, 6), (3, 5), (2, 0)]
captures = []

for move in potential_moves:
    x, y = move
    board[y][x] = 'X'  # Place a black stone
    if is_captured(3, 6, 'O') or is_captured(4, 6, 'O'):
        captures.append(move)
    board[y][x] = '.'  # Reset the board

# Print the result
if captures:
    print(f"Black can capture by playing at: {captures[0]}")
else:
    print("No capture possible.")