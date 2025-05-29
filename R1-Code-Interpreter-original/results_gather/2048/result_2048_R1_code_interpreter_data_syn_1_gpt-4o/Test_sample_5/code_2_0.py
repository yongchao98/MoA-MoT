# Initial matrix
matrix = [
    [8, 32],
    [8, 32]
]

# Move sequence
sequence = 'du'

# Function to perform the moves
def move(matrix, sequence):
    for move in sequence:
        if move == 'u':
            # Move up
            if matrix[1][0] != 0:
                if matrix[0][0] == matrix[1][0]:
                    matrix[0][0] *= 2
                    matrix[1][0] = 0
                elif matrix[0][0] == 0:
                    matrix[0][0], matrix[1][0] = matrix[1][0], 0
        elif move == 'd':
            # Move down
            if matrix[0][0] != 0:
                if matrix[1][0] == matrix[0][0]:
                    matrix[1][0] *= 2
                    matrix[0][0] = 0
                elif matrix[1][0] == 0:
                    matrix[1][0], matrix[0][0] = matrix[0][0], 0
    return matrix

# Perform the moves
result_matrix = move(matrix, sequence)

# Print the result
print(result_matrix)