# Define the board
board = [
    ['.', '.', '.', 'X', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', 'X', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O'],
    ['O', '.', 'X', '.', 'X', 'O', '.', 'X', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Simulate placing a Black stone at E7
board[3][4] = 'X'  # E7 is at row 3, column 4

# Check if the White stone at D7 is captured
def is_captured(board, row, col):
    # Check if the stone at (row, col) is surrounded
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dr, dc in directions:
        r, c = row + dr, col + dc
        if 0 <= r < 10 and 0 <= c < 10 and board[r][c] == '.':
            return False
    return True

# Check if the White stone at D7 is captured
captured = is_captured(board, 3, 3)  # D7 is at row 3, column 3

print("Captured:", captured)