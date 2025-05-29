# Initial matrix
matrix = [
    [0, 0, 0],
    [2, 4, 2],
    [2, 0, 0]
]

# Initial position of the number to move
position = (2, 0)  # (row, column) in zero-indexed

# Move sequence
moves = 'rdl'

# Function to perform the move
def move_number(matrix, position, moves):
    for move in moves:
        row, col = position
        if move == 'r':  # Move right
            if col < 2:  # Check boundary
                if matrix[row][col + 1] == 0:
                    matrix[row][col + 1] = matrix[row][col]
                    matrix[row][col] = 0
                    position = (row, col + 1)
                elif matrix[row][col + 1] == matrix[row][col]:
                    matrix[row][col + 1] *= 2
                    matrix[row][col] = 0
                    position = (row, col + 1)
        elif move == 'd':  # Move down
            if row < 2:  # Check boundary
                if matrix[row + 1][col] == 0:
                    matrix[row + 1][col] = matrix[row][col]
                    matrix[row][col] = 0
                    position = (row + 1, col)
                elif matrix[row + 1][col] == matrix[row][col]:
                    matrix[row + 1][col] *= 2
                    matrix[row][col] = 0
                    position = (row + 1, col)
        elif move == 'l':  # Move left
            if col > 0:  # Check boundary
                if matrix[row][col - 1] == 0:
                    matrix[row][col - 1] = matrix[row][col]
                    matrix[row][col] = 0
                    position = (row, col - 1)
                elif matrix[row][col - 1] == matrix[row][col]:
                    matrix[row][col - 1] *= 2
                    matrix[row][col] = 0
                    position = (row, col - 1)
    return matrix

# Perform the moves
final_matrix = move_number(matrix, position, moves)

# Print the final matrix
print(final_matrix)