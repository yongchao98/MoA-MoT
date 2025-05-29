# Define the board
board = [
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'X', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'O', '.', 'X', 'O', '.', '.', '.'],
    ['.', '.', 'X', 'O', 'X', 'X', 'O', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', 'O', '.', '.'],
    ['.', 'X', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
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
        if cx < 0 or cx >= len(board) or cy < 0 or cy >= len(board[0]) or board[cx][cy] == '.':
            captured = False
            continue

        # If the stone is of the same color, continue exploring
        if board[cx][cy] == color:
            stack.extend([(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)])

    return captured

# Simulate placing a black stone at E5 (index 5, 4)
board[5][4] = 'X'

# Check if the white stone at F5 is captured
captured_stones = 0
if is_captured(board, 5, 5, 'O'):
    captured_stones += 1

print(captured_stones)