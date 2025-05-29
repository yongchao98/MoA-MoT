# Initial matrix
matrix = [
    [32, 8, 32, 0],
    [0, 0, 32, 0],
    [0, 8, 0, 0],
    [4, 0, 0, 0]
]

# Initial position of the number
position = (2, 2)

# Move sequence
moves = 'rudl'

# Function to perform the move
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'r':  # Move right
            while y + 1 < 4 and (matrix[x][y + 1] == 0 or matrix[x][y + 1] == matrix[x][y]):
                if matrix[x][y + 1] == matrix[x][y]:
                    matrix[x][y + 1] *= 2
                    matrix[x][y] = 0
                    break
                matrix[x][y + 1], matrix[x][y] = matrix[x][y], 0
                y += 1
        elif move == 'u':  # Move up
            while x - 1 >= 0 and (matrix[x - 1][y] == 0 or matrix[x - 1][y] == matrix[x][y]):
                if matrix[x - 1][y] == matrix[x][y]:
                    matrix[x - 1][y] *= 2
                    matrix[x][y] = 0
                    break
                matrix[x - 1][y], matrix[x][y] = matrix[x][y], 0
                x -= 1
        elif move == 'd':  # Move down
            while x + 1 < 4 and (matrix[x + 1][y] == 0 or matrix[x + 1][y] == matrix[x][y]):
                if matrix[x + 1][y] == matrix[x][y]:
                    matrix[x + 1][y] *= 2
                    matrix[x][y] = 0
                    break
                matrix[x + 1][y], matrix[x][y] = matrix[x][y], 0
                x += 1
        elif move == 'l':  # Move left
            while y - 1 >= 0 and (matrix[x][y - 1] == 0 or matrix[x][y - 1] == matrix[x][y]):
                if matrix[x][y - 1] == matrix[x][y]:
                    matrix[x][y - 1] *= 2
                    matrix[x][y] = 0
                    break
                matrix[x][y - 1], matrix[x][y] = matrix[x][y], 0
                y -= 1

    return matrix

# Perform the moves
final_matrix = move_number(matrix, position, moves)

# Print the final matrix
print(final_matrix)