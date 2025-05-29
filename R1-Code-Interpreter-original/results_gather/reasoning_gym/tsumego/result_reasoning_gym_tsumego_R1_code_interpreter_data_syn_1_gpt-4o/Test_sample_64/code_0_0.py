# Define the board
board = [
    ['X', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', 'O', 'O', '.', 'O', '.', '.'],
    ['.', '.', '.', 'O', '.', 'X', 'O', 'X', 'O', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['O', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.']
]

# Function to check liberties of a stone
def check_liberties(board, x, y):
    liberties = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(board) and 0 <= ny < len(board[0]) and board[nx][ny] == '.':
            liberties += 1
    return liberties

# Check liberties of the White stones at F8, G7, and G6 after Black plays at G8
board[3][6] = 'X'  # Black plays at G8

liberties_F8 = check_liberties(board, 3, 5)
liberties_G7 = check_liberties(board, 4, 6)
liberties_G6 = check_liberties(board, 5, 6)

print(liberties_F8, liberties_G7, liberties_G6)