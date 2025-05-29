# Initial matrix
matrix = [
    [0, 0, 0, 32],
    [0, 8, 0, 0],
    [0, 0, 16, 4],
    [32, 32, 4, 0]
]

# Initial position of the number to move
row, col = 3, 1  # (4, 2) in 0-indexed is (3, 1)

# Sequence of moves
moves = 'uddd'

# Function to move the number
def move_number(matrix, row, col, moves):
    for move in moves:
        if move == 'u' and row > 0:
            if matrix[row-1][col] == 0:
                matrix[row-1][col] = matrix[row][col]
                matrix[row][col] = 0
            elif matrix[row-1][col] == matrix[row][col]:
                matrix[row-1][col] *= 2
                matrix[row][col] = 0
            row -= 1
        elif move == 'd' and row < 3:
            if matrix[row+1][col] == 0:
                matrix[row+1][col] = matrix[row][col]
                matrix[row][col] = 0
            elif matrix[row+1][col] == matrix[row][col]:
                matrix[row+1][col] *= 2
                matrix[row][col] = 0
            row += 1
        elif move == 'l' and col > 0:
            if matrix[row][col-1] == 0:
                matrix[row][col-1] = matrix[row][col]
                matrix[row][col] = 0
            elif matrix[row][col-1] == matrix[row][col]:
                matrix[row][col-1] *= 2
                matrix[row][col] = 0
            col -= 1
        elif move == 'r' and col < 3:
            if matrix[row][col+1] == 0:
                matrix[row][col+1] = matrix[row][col]
                matrix[row][col] = 0
            elif matrix[row][col+1] == matrix[row][col]:
                matrix[row][col+1] *= 2
                matrix[row][col] = 0
            col += 1
    return matrix

# Apply the moves
result_matrix = move_number(matrix, row, col, moves)

# Print the resulting matrix
print(result_matrix)