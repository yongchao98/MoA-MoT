# Define the board state
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', 'X', 'X', 'X'],
    ['.', '.', '.', '.', '.', 'X', 'X', 'O', 'O'],
    ['.', 'O', '.', '.', '.', 'O', 'X', 'O', '.']
]

# Function to check liberties of a stone
def check_liberties(board, x, y):
    liberties = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 9 and 0 <= ny < 9 and board[nx][ny] == '.':
            liberties += 1
    return liberties

# Check liberties for the critical White stones
liberties_H8 = check_liberties(board, 1, 7)
liberties_I8 = check_liberties(board, 1, 8)
liberties_I2 = check_liberties(board, 7, 8)
liberties_I1 = check_liberties(board, 8, 8)

# Print the liberties
print("Liberties of H8:", liberties_H8)
print("Liberties of I8:", liberties_I8)
print("Liberties of I2:", liberties_I2)
print("Liberties of I1:", liberties_I1)