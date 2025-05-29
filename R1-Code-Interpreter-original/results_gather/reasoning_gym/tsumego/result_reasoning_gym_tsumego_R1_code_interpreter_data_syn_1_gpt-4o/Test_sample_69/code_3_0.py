# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', 'O', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'X', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', 'O', 'O', 'O', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.', '.', '.']
]

# Function to check if a move captures a group
def captures_group(board, move, group):
    x, y = move
    board[y][x] = 'X'  # Simulate the move
    for gx, gy in group:
        # Check if any adjacent position is empty (liberty)
        if any(0 <= gx + dx < len(board[0]) and 0 <= gy + dy < len(board) and board[gy + dy][gx + dx] == '.' 
               for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]):
            return False
    return True

# Define the group at F12, H12
group_F12_H12 = [(5, 1), (7, 1)]

# Check if playing at G12 captures the group
move_G12 = (6, 1)
captures = captures_group(board, move_G12, group_F12_H12)

print(captures)