# Define the board
board = [
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', 'X', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', 'O', 'O', 'X', 'O', '.'],
    ['.', '.', '.', '.', 'X', 'O', '.', '.', 'X', '.'],
    ['O', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', 'O', '.', '.', '.', 'O', '.', '.', '.', 'O'],
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 10 and 0 <= y < 10

# Function to get the liberties of a stone
def get_liberties(x, y, color):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    liberties = 0
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if is_within_board(nx, ny) and board[nx][ny] == '.':
            liberties += 1
    return liberties

# Check liberties of the group at F8 and G8
liberties_F8 = get_liberties(2, 5, 'O')
liberties_G8 = get_liberties(2, 6, 'O')

# Print the total liberties
print(liberties_F8 + liberties_G8)