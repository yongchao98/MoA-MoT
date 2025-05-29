# Initial matrix
matrix = [
    [0, 4, 0],
    [16, 16, 4],
    [0, 32, 0]
]

# Initial position of the number to move
x, y = 1, 1

# Move sequence
moves = 'ull'

# Function to move the number
def move_number(matrix, x, y, moves):
    for move in moves:
        if move == 'u' and x > 0:
            # Move up
            if matrix[x-1][y] == matrix[x][y]:
                matrix[x-1][y] *= 2
                matrix[x][y] = 0
            elif matrix[x-1][y] == 0:
                matrix[x-1][y] = matrix[x][y]
                matrix[x][y] = 0
            x -= 1
        elif move == 'l' and y > 0:
            # Move left
            if matrix[x][y-1] == matrix[x][y]:
                matrix[x][y-1] *= 2
                matrix[x][y] = 0
            elif matrix[x][y-1] == 0:
                matrix[x][y-1] = matrix[x][y]
                matrix[x][y] = 0
            y -= 1
    return matrix

# Apply the moves
result_matrix = move_number(matrix, x, y, moves)

# Print the resulting matrix
print(result_matrix)