# Initial matrix
matrix = [
    [0, 0, 0, 0, 0],
    [8, 2, 4, 0, 0],
    [8, 0, 0, 0, 4],
    [0, 4, 0, 16, 0],
    [0, 0, 0, 0, 4]
]

# Initial position of the number to move
position = (3, 3)

# Movement sequence
sequence = 'lullll'

# Function to move the number according to the sequence
def move_number(matrix, position, sequence):
    x, y = position
    for move in sequence:
        if move == 'l' and y > 0:
            if matrix[x][y-1] == matrix[x][y]:
                matrix[x][y-1] *= 2
                matrix[x][y] = 0
            elif matrix[x][y-1] == 0:
                matrix[x][y-1] = matrix[x][y]
                matrix[x][y] = 0
            y -= 1
        elif move == 'r' and y < 4:
            if matrix[x][y+1] == matrix[x][y]:
                matrix[x][y+1] *= 2
                matrix[x][y] = 0
            elif matrix[x][y+1] == 0:
                matrix[x][y+1] = matrix[x][y]
                matrix[x][y] = 0
            y += 1
        elif move == 'u' and x > 0:
            if matrix[x-1][y] == matrix[x][y]:
                matrix[x-1][y] *= 2
                matrix[x][y] = 0
            elif matrix[x-1][y] == 0:
                matrix[x-1][y] = matrix[x][y]
                matrix[x][y] = 0
            x -= 1
        elif move == 'd' and x < 4:
            if matrix[x+1][y] == matrix[x][y]:
                matrix[x+1][y] *= 2
                matrix[x][y] = 0
            elif matrix[x+1][y] == 0:
                matrix[x+1][y] = matrix[x][y]
                matrix[x][y] = 0
            x += 1
    return matrix

# Apply the movement
final_matrix = move_number(matrix, position, sequence)

# Print the final matrix
print(final_matrix)