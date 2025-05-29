# Initial matrix
matrix = [
    [0, 4, 0, 4],
    [0, 0, 0, 0],
    [0, 4, 16, 0],
    [0, 0, 0, 0]
]

# Position of the number to move
position = (2, 3)  # (row, column) in 0-indexed

# Move sequence
moves = 'duuu'

# Function to move and merge numbers
def move_and_merge(matrix, position, moves):
    row, col = position
    for move in moves:
        if move == 'd' and row < 3:
            if matrix[row + 1][col] == matrix[row][col]:
                matrix[row + 1][col] *= 2
                matrix[row][col] = 0
            elif matrix[row + 1][col] == 0:
                matrix[row + 1][col] = matrix[row][col]
                matrix[row][col] = 0
            row += 1
        elif move == 'u' and row > 0:
            if matrix[row - 1][col] == matrix[row][col]:
                matrix[row - 1][col] *= 2
                matrix[row][col] = 0
            elif matrix[row - 1][col] == 0:
                matrix[row - 1][col] = matrix[row][col]
                matrix[row][col] = 0
            row -= 1
        elif move == 'l' and col > 0:
            if matrix[row][col - 1] == matrix[row][col]:
                matrix[row][col - 1] *= 2
                matrix[row][col] = 0
            elif matrix[row][col - 1] == 0:
                matrix[row][col - 1] = matrix[row][col]
                matrix[row][col] = 0
            col -= 1
        elif move == 'r' and col < 3:
            if matrix[row][col + 1] == matrix[row][col]:
                matrix[row][col + 1] *= 2
                matrix[row][col] = 0
            elif matrix[row][col + 1] == 0:
                matrix[row][col + 1] = matrix[row][col]
                matrix[row][col] = 0
            col += 1
    return matrix

# Apply the moves
result_matrix = move_and_merge(matrix, position, moves)

# Print the result
print(result_matrix)