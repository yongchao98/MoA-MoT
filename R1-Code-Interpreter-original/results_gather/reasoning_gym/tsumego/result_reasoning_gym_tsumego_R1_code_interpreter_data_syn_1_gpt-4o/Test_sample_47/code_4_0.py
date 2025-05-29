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

# Simulate placing a Black stone at G2
board[8][6] = 'X'  # G2 is at row 8, column 6

# Check if the White stones at F3, G3, and H3 are captured
def is_captured(board, row, col):
    # Check if the stone at (row, col) is surrounded
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dr, dc in directions:
        r, c = row + dr, col + dc
        if 0 <= r < 10 and 0 <= c < 10 and board[r][c] == '.':
            return False
    return True

# Check if the White stones at F3, G3, and H3 are captured
captured_F3 = is_captured(board, 7, 5)  # F3 is at row 7, column 5
captured_G3 = is_captured(board, 7, 6)  # G3 is at row 7, column 6
captured_H3 = is_captured(board, 7, 7)  # H3 is at row 7, column 7

print("Captured F3:", captured_F3)
print("Captured G3:", captured_G3)
print("Captured H3:", captured_H3)