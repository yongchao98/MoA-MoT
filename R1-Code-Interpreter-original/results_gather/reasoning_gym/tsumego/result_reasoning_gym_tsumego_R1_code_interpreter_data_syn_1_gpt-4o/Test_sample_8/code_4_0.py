# Define the board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', 'O', '.'],
    ['.', 'X', 'O', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 12 and 0 <= y < 12

# Function to get the liberties of a stone
def get_liberties(x, y):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    liberties = 0
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if is_within_board(nx, ny) and board[nx][ny] == '.':
            liberties += 1
    return liberties

# Check the liberties of the White stone at G11 (row 11, column 7)
x, y = 1, 6  # G11
liberties = get_liberties(x, y)
print("Liberties of the White stone at G11:", liberties)