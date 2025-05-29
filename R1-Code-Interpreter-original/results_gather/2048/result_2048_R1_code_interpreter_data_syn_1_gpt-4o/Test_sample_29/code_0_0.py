# Initial matrix
matrix = [
    [16, 4, 0],
    [0, 2, 0],
    [32, 2, 4]
]

# Initial position of the number to move
row, col = 2, 2  # 0-indexed, corresponds to (2, 3) in 1-indexed

# Sequence of moves
moves = 'lrr'

# Function to perform the move
def move_number(matrix, row, col, direction):
    if direction == 'l' and col > 0:
        # Move left
        if matrix[row][col - 1] == matrix[row][col]:
            matrix[row][col - 1] *= 2
            matrix[row][col] = 0
        elif matrix[row][col - 1] == 0:
            matrix[row][col - 1] = matrix[row][col]
            matrix[row][col] = 0
        col -= 1
    elif direction == 'r' and col < 2:
        # Move right
        if matrix[row][col + 1] == matrix[row][col]:
            matrix[row][col + 1] *= 2
            matrix[row][col] = 0
        elif matrix[row][col + 1] == 0:
            matrix[row][col + 1] = matrix[row][col]
            matrix[row][col] = 0
        col += 1
    return row, col

# Apply the sequence of moves
for move in moves:
    row, col = move_number(matrix, row, col, move)

# Print the final matrix
print(matrix)