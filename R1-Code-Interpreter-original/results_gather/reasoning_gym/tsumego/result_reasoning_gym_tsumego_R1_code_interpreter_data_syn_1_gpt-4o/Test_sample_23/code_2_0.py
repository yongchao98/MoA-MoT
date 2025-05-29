# Define the board
board = [
    ['.', '.', '.', '.', 'X', '.', '.', 'X', 'X'],
    ['.', '.', '.', '.', '.', 'X', 'X', 'O', 'O'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', '.'],
    ['.', 'O', 'O', '.', '.', '.', 'X', 'X', 'O'],
    ['.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.']
]

# Function to check liberties of a stone
def check_liberties(x, y, board):
    liberties = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 9 and 0 <= ny < 9 and board[ny][nx] == '.':
            liberties += 1
    return liberties

# Check liberties for H8 and H6
liberties_H8 = check_liberties(7, 7, board)
liberties_H6 = check_liberties(7, 5, board)

# Print the liberties
print("Liberties of H8:", liberties_H8)
print("Liberties of H6:", liberties_H6)