# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', 'X', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'X', '.', 'X', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', 'X', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', 'O', '.', '.', '.', '.']
]

# Function to check liberties of a stone
def check_liberties(board, x, y):
    liberties = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 12 and 0 <= ny < 12 and board[nx][ny] == '.':
            liberties += 1
    return liberties

# Place a Black stone at D6 (index 6, 3)
board[6][3] = 'X'

# Check liberties of White stones at D5, D4, D3
liberties_D5 = check_liberties(board, 7, 3)
liberties_D4 = check_liberties(board, 8, 3)
liberties_D3 = check_liberties(board, 9, 3)

# Print the liberties
print(liberties_D5, liberties_D4, liberties_D3)