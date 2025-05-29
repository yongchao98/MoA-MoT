# Simulate the move at D2 and check the liberties of the group at C3, D3, C2
def check_liberties(board, x, y):
    # Check the liberties of the stone at (x, y)
    liberties = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(board) and 0 <= ny < len(board[0]) and board[nx][ny] == '.':
            liberties += 1
    return liberties

# Board representation
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'X', 'X', 'X', '.', '.', 'X', 'X', '.'],
    ['.', 'X', 'O', 'O', 'X', 'O', 'X', '.', 'O'],
    ['.', 'X', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'O', 'X', '.', '.', '.', '.', '.']
]

# Simulate the move at D2
board[7][3] = 'X'  # D2 is (7, 3) in 0-indexed

# Check liberties of the group at C3, D3, C2
liberties_C3 = check_liberties(board, 6, 2)
liberties_D3 = check_liberties(board, 6, 3)
liberties_C2 = check_liberties(board, 7, 2)

# Output the total liberties
total_liberties = liberties_C3 + liberties_D3 + liberties_C2
print(total_liberties)