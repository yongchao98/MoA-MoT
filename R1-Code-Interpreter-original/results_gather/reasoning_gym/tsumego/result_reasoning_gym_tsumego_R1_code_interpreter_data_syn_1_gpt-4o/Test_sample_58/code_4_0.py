# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', '.', '.', '.', '.', 'O', '.', '.'],
    ['O', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['O', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'X', 'X', 'X', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.']
]

# Function to check if a group of stones is captured
def is_captured(board, x, y, color):
    visited = set()
    stack = [(x, y)]
    captured = True

    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))

        # Check if the stone is on the edge or has a liberty
        if cx < 0 or cx >= 10 or cy < 0 or cy >= 10 or board[cy][cx] == '.':
            return False

        # Check adjacent stones
        if board[cy][cx] == color:
            stack.extend([(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)])

    return captured

# Simulate the move at C5
board[5][2] = 'X'

# Check if the group is captured
captured_stones = 0
for y in range(10):
    for x in range(10):
        if board[y][x] == 'O' and is_captured(board, x, y, 'O'):
            captured_stones += 1

print(captured_stones)