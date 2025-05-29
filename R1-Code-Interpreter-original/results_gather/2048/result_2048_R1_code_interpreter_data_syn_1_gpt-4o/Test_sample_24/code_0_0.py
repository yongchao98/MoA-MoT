# Initial matrix
matrix = [
    [0, 0, 0],
    [2, 4, 2],
    [2, 0, 0]
]

# Initial position of the number to move
position = (2, 0)  # (row, column) in 0-indexed

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
        elif move == 'd':  # Move down
            if row < 2:  # Check boundary
                if matrix[row + 1][col] == 0:
                    matrix[row + 1][col] = matrix[row][col]
                    matrix[row][col] = 0
                    position = (row + 1, col)
        elif move == 'l':  # Move left
            if col > 0:  # Check boundary
                if matrix[row][col - 1] == 0:
                    matrix[row][col - 1] = matrix[row][col]
                    matrix[row][col] = 0
                    position = (row, col - 1)
        elif move == 'u':  # Move up
            if row > 0:  # Check boundary
                if matrix[row - 1][col] == 0:
                    matrix[row - 1][col] = matrix[row][col]
                    matrix[row][col] = 0
                    position = (row - 1, col)
    return matrix

# Perform the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)