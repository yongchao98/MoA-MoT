# Initial matrix
matrix = [
    [16, 8, 0, 2],
    [0, 0, 32, 0],
    [0, 4, 32, 0],
    [0, 32, 0, 2]
]

# Initial position of the number to move
position = (3, 3)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'ddud'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    for move in moves:
        row, col = position
        if move == 'u' and row > 0:
            if matrix[row-1][col] == 0:
                matrix[row-1][col] = matrix[row][col]
                matrix[row][col] = 0
                position = (row-1, col)
            elif matrix[row-1][col] == matrix[row][col]:
                matrix[row-1][col] *= 2
                matrix[row][col] = 0
                position = (row-1, col)
        elif move == 'd' and row < 3:
            if matrix[row+1][col] == 0:
                matrix[row+1][col] = matrix[row][col]
                matrix[row][col] = 0
                position = (row+1, col)
            elif matrix[row+1][col] == matrix[row][col]:
                matrix[row+1][col] *= 2
                matrix[row][col] = 0
                position = (row+1, col)
        elif move == 'l' and col > 0:
            if matrix[row][col-1] == 0:
                matrix[row][col-1] = matrix[row][col]
                matrix[row][col] = 0
                position = (row, col-1)
            elif matrix[row][col-1] == matrix[row][col]:
                matrix[row][col-1] *= 2
                matrix[row][col] = 0
                position = (row, col-1)
        elif move == 'r' and col < 3:
            if matrix[row][col+1] == 0:
                matrix[row][col+1] = matrix[row][col]
                matrix[row][col] = 0
                position = (row, col+1)
            elif matrix[row][col+1] == matrix[row][col]:
                matrix[row][col+1] *= 2
                matrix[row][col] = 0
                position = (row, col+1)
    return matrix

# Apply the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)