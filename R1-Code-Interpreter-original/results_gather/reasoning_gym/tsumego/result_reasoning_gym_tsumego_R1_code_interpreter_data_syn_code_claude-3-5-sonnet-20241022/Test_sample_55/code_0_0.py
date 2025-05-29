def count_liberties(board, x, y):
    liberties = 0
    # Check adjacent positions (up, right, down, left)
    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        new_x, new_y = x + dx, y + dy
        if 0 <= new_x < 10 and 0 <= new_y < 10:
            if board[new_y][new_x] == '.':
                liberties += 1

# Create the board
board = [['.'] * 10 for _ in range(10)]
# Fill in the stones from the given position
position = """
. . . . . O . . . .
. . . X . X X . . .
. . . . X O O X O .
. . . . X O . . X .
O . . . . X . . . .
. . . . . . . . O .
. O . . . O . . . O
. . . O . . . . . .
. . . . . . . . . O
. . . . . . . . . .
"""

# Convert the position string to board state
rows = position.strip().split('\n')
for i, row in enumerate(rows):
    stones = row.split()
    for j, stone in enumerate(stones):
        board[i][j] = stone

# Check liberties of white stones at F8 and G8
f8_liberties = count_liberties(board, 5, 2)  # F8
g8_liberties = count_liberties(board, 6, 2)  # G8

# Check if G7 is the capturing move
g7_x, g7_y = 6, 3  # G7 coordinates
is_empty = board[g7_y][g7_x] == '.'

print(f"F8 liberties: {f8_liberties}")
print(f"G8 liberties: {g8_liberties}")
print(f"G7 is empty: {is_empty}")
print(f"G7 is capturing move: {is_empty and f8_liberties == 1 and g8_liberties == 1}")