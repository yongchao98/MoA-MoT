def count_liberties(board, x, y):
    liberties = 0
    # Check adjacent positions (up, right, down, left)
    for dx, dy in [(0, -1), (1, 0), (0, 1), (-1, 0)]:  # Fixed directions
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 10 and 0 <= new_y < 10:
            if board[new_y][new_x] == '.':
                liberties += 1
    return liberties

# Create the board (10x10)
board = [
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],  # 10
    ['.', '.', '.', 'X', '.', 'X', 'X', '.', '.', '.'],  # 9
    ['.', '.', '.', '.', 'X', 'O', 'O', 'X', 'O', '.'],  # 8
    ['.', '.', '.', '.', 'X', 'O', '.', '.', 'X', '.'],  # 7
    ['O', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],  # 6
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],  # 5
    ['.', 'O', '.', '.', '.', 'O', '.', '.', '.', 'O'],  # 4
    ['.', '.', '.', 'O', '.', '.', '.', '.', '.', '.'],  # 3
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'O'],  # 2
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],  # 1
]

# Check liberties of white stones at F8 and G8
f8_x, f8_y = 5, 2  # F8 coordinates (0-based)
g8_x, g8_y = 6, 2  # G8 coordinates (0-based)
g7_x, g7_y = 6, 3  # G7 coordinates (0-based)

f8_liberties = count_liberties(board, f8_x, f8_y)
g8_liberties = count_liberties(board, g8_x, g8_y)
is_empty = board[g7_y][g7_x] == '.'

print(f"F8 liberties: {f8_liberties}")
print(f"G8 liberties: {g8_liberties}")
print(f"G7 is empty: {is_empty}")
print(f"Connected white group at F8-G8 can be captured by G7: {is_empty}")

# Additional check for the white stone at F7
f7_x, f7_y = 5, 3  # F7 coordinates
f7_liberties = count_liberties(board, f7_x, f7_y)
print(f"F7 liberties: {f7_liberties}")