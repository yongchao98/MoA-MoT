# Define the board as a 2D list
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', 'O', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check liberties of a stone at a given position
def check_liberties(board, x, y):
    liberties = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(board) and 0 <= ny < len(board[0]) and board[nx][ny] == '.':
            liberties += 1
    return liberties

# Play at H10
board[2][7] = 'X'  # H10 is at (2, 7) in 0-indexed

# Check liberties of the White stone at G9 after playing H10
liberties_G9 = check_liberties(board, 3, 6)  # G9 is at (3, 6) in 0-indexed

print("Liberties of White stone at G9 after playing H10:", liberties_G9)