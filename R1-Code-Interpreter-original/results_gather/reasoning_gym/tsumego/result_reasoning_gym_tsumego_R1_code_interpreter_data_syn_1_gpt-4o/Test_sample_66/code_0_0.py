# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
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

# Check the move at F11
x, y = 0, 5  # F11
if board[x][y] == '.':
    # Check the liberties of the white stone at E11
    liberties = check_liberties(x, y)
    print("Liberties for E11 after F11 move:", liberties)

# Check if the move captures any stones
captured = False
if liberties == 0:
    captured = True

print("Move F11 captures stones:", captured)