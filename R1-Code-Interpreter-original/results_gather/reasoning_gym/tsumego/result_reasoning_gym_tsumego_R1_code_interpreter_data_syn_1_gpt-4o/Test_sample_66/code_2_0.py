# Redefine the board
board = [
    ['.', '.', '.', '.', '.', '.', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'X'],
    ['.', '.', '.', 'X', '.', '.', '.', 'X', 'O', 'O', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', 'X', 'O', 'X'],
    ['X', '.', '.', '.', '.', '.', '.', '.', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', 'X', '.', '.']
]

# Function to check liberties
def check_liberties(x, y):
    liberties = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 11 and 0 <= ny < 11 and board[nx][ny] == '.':
            liberties += 1
    return liberties

# Check the move at H6
x, y = 5, 7  # H6
if board[x][y] == '.':
    # Check the liberties of the white group at G5, H5, and H4
    liberties_g5 = check_liberties(6, 6)
    liberties_h5 = check_liberties(6, 7)
    liberties_h4 = check_liberties(7, 7)
    print("Liberties for G5 after H6 move:", liberties_g5)
    print("Liberties for H5 after H6 move:", liberties_h5)
    print("Liberties for H4 after H6 move:", liberties_h4)

# Check if the move captures any stones
captured_g5 = liberties_g5 == 0
captured_h5 = liberties_h5 == 0
captured_h4 = liberties_h4 == 0

print("Move H6 captures stones at G5:", captured_g5)
print("Move H6 captures stones at H5:", captured_h5)
print("Move H6 captures stones at H4:", captured_h4)