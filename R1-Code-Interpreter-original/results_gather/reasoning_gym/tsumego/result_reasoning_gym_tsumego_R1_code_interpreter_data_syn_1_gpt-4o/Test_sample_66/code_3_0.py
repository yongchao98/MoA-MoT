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

# Check the move at J2
x, y = 9, 8  # J2
if board[x][y] == '.':
    # Check the liberties of the white stone at J3
    liberties_j3 = check_liberties(8, 9)
    print("Liberties for J3 after J2 move:", liberties_j3)

# Check if the move captures any stones
captured_j3 = liberties_j3 == 0

print("Move J2 captures stones at J3:", captured_j3)