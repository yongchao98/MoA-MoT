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

# Function to check if a specific group is captured
def is_group_captured(board, group):
    for x, y in group:
        # Check all adjacent positions for liberties
        for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            nx, ny = x + dx, y + dy
            if 0 <= nx < 10 and 0 <= ny < 10 and board[ny][nx] == '.':
                return False
    return True

# Define the group at B9, C9, D9
group = [(1, 1), (2, 1), (3, 1)]

# Check if playing at B8 captures this group
move = (1, 2)  # B8
board[move[1]][move[0]] = 'X'  # Simulate the move
captured = is_group_captured(board, group)
captured_stones = len(group) if captured else 0

print(captured_stones)