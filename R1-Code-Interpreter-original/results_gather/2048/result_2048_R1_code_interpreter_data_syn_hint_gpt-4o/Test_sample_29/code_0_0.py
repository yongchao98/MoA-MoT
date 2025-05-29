# Initial matrix
matrix = [
    [16, 4, 0],
    [0, 2, 0],
    [32, 2, 4]
]

# Initial position of the number to move
row, col = 2, 3

# Sequence of moves
moves = 'lrr'

# Function to perform the move
def move_number(matrix, row, col, moves):
    for move in moves:
        if move == 'l' and col > 1:
            # Move left
            if matrix[row-1][col-2] == matrix[row-1][col-1]:
                matrix[row-1][col-2] *= 2
                matrix[row-1][col-1] = 0
            else:
                matrix[row-1][col-2], matrix[row-1][col-1] = matrix[row-1][col-1], matrix[row-1][col-2]
            col -= 1
        elif move == 'r' and col < 3:
            # Move right
            if matrix[row-1][col] == matrix[row-1][col-1]:
                matrix[row-1][col] *= 2
                matrix[row-1][col-1] = 0
            else:
                matrix[row-1][col], matrix[row-1][col-1] = matrix[row-1][col-1], matrix[row-1][col]
            col += 1
    return matrix

# Perform the moves
result_matrix = move_number(matrix, row, col, moves)

# Print the resulting matrix
print(result_matrix)