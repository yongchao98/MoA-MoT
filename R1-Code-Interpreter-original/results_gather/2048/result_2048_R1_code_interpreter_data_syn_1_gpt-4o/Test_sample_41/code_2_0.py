# Initial matrix
matrix = [
    [0, 4, 16, 0],
    [0, 0, 2, 0],
    [4, 0, 0, 0],
    [2, 0, 0, 0]
]

# Initial position of the number to move
position = (3, 0)  # (row, column) in 0-indexed

# Sequence of moves
moves = 'duuu'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    for move in moves:
        row, col = position
        if move == 'd' and row < 3:
            # Move down
            if matrix[row + 1][col] == 0:
                matrix[row + 1][col] = matrix[row][col]
                matrix[row][col] = 0
            position = (row + 1, col)
        elif move == 'u' and row > 0:
            # Move up
            if matrix[row - 1][col] == 0:
                matrix[row - 1][col] = matrix[row][col]
                matrix[row][col] = 0
            position = (row - 1, col)
        elif move == 'l' and col > 0:
            # Move left
            if matrix[row][col - 1] == 0:
                matrix[row][col - 1] = matrix[row][col]
                matrix[row][col] = 0
            position = (row, col - 1)
        elif move == 'r' and col < 3:
            # Move right
            if matrix[row][col + 1] == 0:
                matrix[row][col + 1] = matrix[row][col]
                matrix[row][col] = 0
            position = (row, col + 1)
    return matrix

# Apply the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)