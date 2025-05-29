# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', 'O', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', 'X', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'O', 'O', 'O', 'O', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', '.', 'X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', '.', '.', '.', '.', 'X', '.', 'X', '.']
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

# Check liberties of the White stone at E8 after placing a Black stone at E7
board[5][4] = 'X'  # Place Black stone at E7
liberties_E8 = check_liberties(board, 5, 4)

# Output the liberties of the White stone at E8
print(liberties_E8)