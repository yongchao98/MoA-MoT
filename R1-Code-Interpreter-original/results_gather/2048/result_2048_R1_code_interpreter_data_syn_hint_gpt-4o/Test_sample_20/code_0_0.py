# Initial matrix
matrix = [
    [0, 2, 0],
    [0, 32, 0],
    [0, 8, 0]
]

# Initial position of the number to move
position = (2, 2)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'ddu'

# Function to perform a move
def move_number(matrix, position, direction):
    row, col = position
    if direction == 'd' and row < 2:
        # Move down
        if matrix[row + 1][col] == matrix[row][col]:
            matrix[row + 1][col] *= 2
            matrix[row][col] = 0
        elif matrix[row + 1][col] == 0:
            matrix[row + 1][col] = matrix[row][col]
            matrix[row][col] = 0
        row += 1
    elif direction == 'u' and row > 0:
        # Move up
        if matrix[row - 1][col] == matrix[row][col]:
            matrix[row - 1][col] *= 2
            matrix[row][col] = 0
        elif matrix[row - 1][col] == 0:
            matrix[row - 1][col] = matrix[row][col]
            matrix[row][col] = 0
        row -= 1
    return (row, col)

# Perform the sequence of moves
for move in moves:
    position = move_number(matrix, position, move)

# Print the final matrix
print(matrix)