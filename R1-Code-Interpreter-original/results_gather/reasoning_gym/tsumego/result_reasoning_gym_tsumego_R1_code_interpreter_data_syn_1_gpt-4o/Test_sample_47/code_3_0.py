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

# Simulate placing a Black stone at G4
board[6][3] = 'X'  # G4 is at row 6, column 3

# Check if the White stones at H4 and I4 are captured
def is_captured(board, row, col):
    # Check if the stone at (row, col) is surrounded
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dr, dc in directions:
        r, c = row + dr, col + dc
        if 0 <= r < 10 and 0 <= c < 10 and board[r][c] == '.':
            return False
    return True

# Check if the White stones at H4 and I4 are captured
captured_H4 = is_captured(board, 6, 7)  # H4 is at row 6, column 7
captured_I4 = is_captured(board, 6, 8)  # I4 is at row 6, column 8

print("Captured H4:", captured_H4)
print("Captured I4:", captured_I4)