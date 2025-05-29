# Define the board
board = [
    ['.', '.', 'X', '.', '.', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['O', '.', 'X', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'X', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'O', 'O', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', 'X', 'O', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Check the liberties of the White stones at C3, D3, and C2
def check_liberties(board, row, col):
    liberties = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dr, dc in directions:
        r, c = row + dr, col + dc
        if 0 <= r < len(board) and 0 <= c < len(board[0]) and board[r][c] == '.':
            liberties += 1
    return liberties

# Check the liberties of the White stones
liberties_C3 = check_liberties(board, 10, 2)
liberties_D3 = check_liberties(board, 10, 3)
liberties_C2 = check_liberties(board, 11, 2)

# Check the effect of placing a Black stone at D2
board[11][3] = 'X'
liberties_C3_after = check_liberties(board, 10, 2)
liberties_D3_after = check_liberties(board, 10, 3)
liberties_C2_after = check_liberties(board, 11, 2)

print("Liberties before move:")
print(f"C3: {liberties_C3}, D3: {liberties_D3}, C2: {liberties_C2}")
print("Liberties after move at D2:")
print(f"C3: {liberties_C3_after}, D3: {liberties_D3_after}, C2: {liberties_C2_after}")