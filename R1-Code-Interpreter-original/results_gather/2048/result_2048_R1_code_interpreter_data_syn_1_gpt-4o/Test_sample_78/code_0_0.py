# Initial matrix
matrix = [
    [32, 0, 0, 0, 0],
    [0, 4, 4, 2, 0],
    [0, 2, 0, 0, 0],
    [0, 8, 0, 0, 0],
    [0, 0, 0, 0, 0]
]

# Initial position of the number to move
position = (4, 1)  # (row, column) in 0-indexed

# Move sequence
moves = 'lrdlr'

# Function to perform a move
def move(position, direction):
    row, col = position
    if direction == 'l' and col > 0:
        col -= 1
    elif direction == 'r' and col < 4:
        col += 1
    elif direction == 'u' and row > 0:
        row -= 1
    elif direction == 'd' and row < 4:
        row += 1
    return (row, col)

# Simulate the moves
for move_direction in moves:
    position = move(position, move_direction)

# Print the final matrix
print(matrix)