# Initial matrix
matrix = [
    [0, 0, 0, 32],
    [0, 8, 0, 0],
    [0, 0, 16, 4],
    [32, 32, 4, 0]
]

# Initial position of the number to move
position = (3, 1)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'uddd'

# Function to perform a move
def move_number(matrix, position, direction):
    row, col = position
    if direction == 'u' and row > 0:
        target_row = row - 1
        if matrix[target_row][col] == matrix[row][col]:
            matrix[target_row][col] *= 2
            matrix[row][col] = 0
        elif matrix[target_row][col] == 0:
            matrix[target_row][col] = matrix[row][col]
            matrix[row][col] = 0
        position = (target_row, col)
    elif direction == 'd' and row < 3:
        target_row = row + 1
        if matrix[target_row][col] == matrix[row][col]:
            matrix[target_row][col] *= 2
            matrix[row][col] = 0
        elif matrix[target_row][col] == 0:
            matrix[target_row][col] = matrix[row][col]
            matrix[row][col] = 0
        position = (target_row, col)
    return position

# Perform the sequence of moves
for move in moves:
    position = move_number(matrix, position, move)

# Print the final matrix
print(matrix)