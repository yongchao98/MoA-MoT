# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['O', '.', 'X', 'O', 'X', '.', '.', '.', 'X'],
    ['O', 'X', 'O', 'O', '.', '.', '.', '.', '.'],
    ['.', 'X', 'X', 'O', 'X', 'X', '.', '.', '.'],
    ['X', 'X', '.', 'X', '.', 'O', '.', '.', '.']
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

# Simulate the move E3 for Black
board[6][4] = 'X'  # E3

# Check liberties of White stones at D3, C3, and D4
liberties_D3 = check_liberties(board, 6, 3)
liberties_C3 = check_liberties(board, 6, 2)
liberties_D4 = check_liberties(board, 5, 3)

print("Liberties after move E3:")
print("D3:", liberties_D3)
print("C3:", liberties_C3)
print("D4:", liberties_D4)