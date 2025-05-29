# Define the board
board = [
    ['.', 'X', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', 'O', 'X', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'X', '.', 'O', '.', 'X', '.', '.'],
    ['X', '.', 'X', '.', '.', '.', '.', '.', 'X', '.'],
    ['X', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.']
]

# Function to check if a move captures a specific group
def captures_group(board, move, group):
    x, y = move
    # Simulate the move
    board[y][x] = 'X'
    # Check if the group has any liberties left
    for gx, gy in group:
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = gx + dx, gy + dy
            if 0 <= nx < 10 and 0 <= ny < 10 and board[ny][nx] == '.':
                return False
    return True

# Define the group at B9, C9, D9
group = [(1, 1), (2, 1), (3, 1)]

# Check if playing at B8 captures this group
move = (1, 2)  # B8
captured = captures_group(board, move, group)
captured_stones = len(group) if captured else 0

print(captured_stones)