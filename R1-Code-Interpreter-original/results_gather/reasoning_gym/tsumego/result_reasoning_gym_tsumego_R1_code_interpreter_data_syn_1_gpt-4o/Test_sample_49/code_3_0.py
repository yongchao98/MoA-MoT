# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', 'O', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', 'X', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'O', 'O', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', 'X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', 'X', '.', 'X', '.']
]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 13 and 0 <= y < 13

# Function to get the adjacent positions
def get_adjacent(x, y):
    return [(x-1, y), (x+1, y), (x, y-1), (x, y+1)]

# Function to check if a group of stones is captured
def is_captured(x, y, color):
    visited = set()
    stack = [(x, y)]
    captured = True

    while stack:
        cx, cy = stack.pop()
        if (cx, cy) in visited:
            continue
        visited.add((cx, cy))

        for nx, ny in get_adjacent(cx, cy):
            if is_within_board(nx, ny):
                if board[nx][ny] == '.':
                    captured = False
                elif board[nx][ny] == color:
                    stack.append((nx, ny))

    return captured

# Simulate placing a black stone at E9 (index 9, 4)
board[9][4] = 'X'

# Check if any white stones are captured
captured_stones = 0
for x in range(13):
    for y in range(13):
        if board[x][y] == 'O' and is_captured(x, y, 'O'):
            captured_stones += 1

print(captured_stones)