# Initial matrix
matrix = [
    [0, 32, 0, 0],
    [4, 4, 0, 2],
    [16, 0, 0, 32],
    [0, 16, 0, 0]
]

# Initial position of the number to move
row, col = 1, 2  # (2, 3) in 0-indexed

# Sequence of moves
moves = 'ddll'

# Function to perform a move
def move_number(matrix, row, col, direction):
    if direction == 'd' and row < 3:
        # Move down
        if matrix[row + 1][col] == matrix[row][col]:
            matrix[row + 1][col] *= 2
            matrix[row][col] = 0
        elif matrix[row + 1][col] == 0:
            matrix[row + 1][col] = matrix[row][col]
            matrix[row][col] = 0
        row += 1
    elif direction == 'l' and col > 0:
        # Move left
        if matrix[row][col - 1] == matrix[row][col]:
            matrix[row][col - 1] *= 2
            matrix[row][col] = 0
        elif matrix[row][col - 1] == 0:
            matrix[row][col - 1] = matrix[row][col]
            matrix[row][col] = 0
        col -= 1
    return row, col

# Perform the moves
for move in moves:
    row, col = move_number(matrix, row, col, move)

# Print the final matrix
print(matrix)