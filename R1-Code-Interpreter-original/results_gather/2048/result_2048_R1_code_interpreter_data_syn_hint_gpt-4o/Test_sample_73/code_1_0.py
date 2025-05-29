# Initial matrix
matrix = [
    [0, 0, 8, 0, 0],
    [0, 0, 32, 0, 32],
    [0, 32, 0, 0, 4],
    [0, 0, 0, 0, 4],
    [0, 0, 0, 32, 0]
]

# Initial position of the number
position = (4, 4)  # Adjusted to zero-based index

# Sequence of moves
moves = 'lddud'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'l' and y > 0:
            if matrix[x][y-1] == 0:
                matrix[x][y-1] = matrix[x][y]
                matrix[x][y] = 0
            elif matrix[x][y-1] == matrix[x][y]:
                matrix[x][y-1] *= 2
                matrix[x][y] = 0
            y -= 1
        elif move == 'r' and y < 4:
            if matrix[x][y+1] == 0:
                matrix[x][y+1] = matrix[x][y]
                matrix[x][y] = 0
            elif matrix[x][y+1] == matrix[x][y]:
                matrix[x][y+1] *= 2
                matrix[x][y] = 0
            y += 1
        elif move == 'u' and x > 0:
            if matrix[x-1][y] == 0:
                matrix[x-1][y] = matrix[x][y]
                matrix[x][y] = 0
            elif matrix[x-1][y] == matrix[x][y]:
                matrix[x-1][y] *= 2
                matrix[x][y] = 0
            x -= 1
        elif move == 'd' and x < 4:
            if matrix[x+1][y] == 0:
                matrix[x+1][y] = matrix[x][y]
                matrix[x][y] = 0
            elif matrix[x+1][y] == matrix[x][y]:
                matrix[x+1][y] *= 2
                matrix[x][y] = 0
            x += 1
    return matrix

# Update the matrix
updated_matrix = move_number(matrix, position, moves)

# Print the updated matrix
print(updated_matrix)