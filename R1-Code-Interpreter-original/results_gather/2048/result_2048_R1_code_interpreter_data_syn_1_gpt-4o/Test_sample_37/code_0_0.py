# Initial matrix
matrix = [
    [0, 0, 0],
    [4, 32, 0],
    [0, 4, 0]
]

# Initial position of the number to move
start_position = (2, 1)  # (row, column) in 0-indexed

# Move sequence
moves = 'uuu'

# Function to move the number in the matrix
def move_number(matrix, start_position, moves):
    row, col = start_position
    for move in moves:
        if move == 'u' and row > 0:
            # Check if the number can combine with the number above
            if matrix[row-1][col] == matrix[row][col]:
                matrix[row-1][col] *= 2
                matrix[row][col] = 0
            elif matrix[row-1][col] == 0:
                matrix[row-1][col] = matrix[row][col]
                matrix[row][col] = 0
            row -= 1
    return matrix

# Apply the moves
final_matrix = move_number(matrix, start_position, moves)

# Print the final matrix
print(final_matrix)