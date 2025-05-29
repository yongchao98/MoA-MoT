# Initial matrix
matrix = [
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 16],
    [32, 8, 0, 0, 0],
    [4, 0, 0, 16, 0],
    [8, 0, 0, 0, 0]
]

# Initial position of the number to move
start_row, start_col = 4, 3

# Sequence of moves
moves = 'ddddd'

# Simulate the moves
current_row, current_col = start_row, start_col
for move in moves:
    if move == 'd' and current_row < 4:
        current_row += 1

# The number to move is 0, so the matrix remains unchanged
# Print the final matrix
for row in matrix:
    print(row)