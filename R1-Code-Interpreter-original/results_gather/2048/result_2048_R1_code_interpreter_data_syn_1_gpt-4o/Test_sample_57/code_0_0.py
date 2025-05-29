# Initial matrix
matrix = [
    [2, 2, 0, 0],
    [0, 16, 0, 0],
    [4, 8, 2, 0],
    [0, 32, 0, 0]
]

# Initial position of the number to move
row, col = 3, 2

# Movement sequence
moves = 'dull'

# Function to perform the moves
def move_number(matrix, row, col, moves):
    for move in moves:
        if move == 'd':  # Move down
            if row < 3 and matrix[row + 1][col] == 0:
                matrix[row + 1][col] = matrix[row][col]
                matrix[row][col] = 0
                row += 1
        elif move == 'u':  # Move up
            if row > 0 and matrix[row - 1][col] == 0:
                matrix[row - 1][col] = matrix[row][col]
                matrix[row][col] = 0
                row -= 1
        elif move == 'l':  # Move left
            if col > 0 and matrix[row][col - 1] == 0:
                matrix[row][col - 1] = matrix[row][col]
                matrix[row][col] = 0
                col -= 1
        elif move == 'r':  # Move right
            if col < 3 and matrix[row][col + 1] == 0:
                matrix[row][col + 1] = matrix[row][col]
                matrix[row][col] = 0
                col += 1
    return matrix

# Perform the moves
final_matrix = move_number(matrix, row, col, moves)

# Print the final matrix
print(final_matrix)