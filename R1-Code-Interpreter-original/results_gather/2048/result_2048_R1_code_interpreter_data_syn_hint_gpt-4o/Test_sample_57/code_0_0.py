# Initial matrix
matrix = [
    [2, 2, 0, 0],
    [0, 16, 0, 0],
    [4, 8, 2, 0],
    [0, 32, 0, 0]
]

# Initial position of the number to move
position = (3, 2)

# Movement sequence
sequence = 'dull'

# Function to move the number according to the sequence
def move_number(matrix, position, sequence):
    x, y = position
    for move in sequence:
        if move == 'd':  # Move down
            if x < 3:
                if matrix[x+1][y] == matrix[x][y]:
                    matrix[x+1][y] *= 2
                    matrix[x][y] = 0
                elif matrix[x+1][y] == 0:
                    matrix[x+1][y] = matrix[x][y]
                    matrix[x][y] = 0
                x += 1
        elif move == 'u':  # Move up
            if x > 0:
                if matrix[x-1][y] == matrix[x][y]:
                    matrix[x-1][y] *= 2
                    matrix[x][y] = 0
                elif matrix[x-1][y] == 0:
                    matrix[x-1][y] = matrix[x][y]
                    matrix[x][y] = 0
                x -= 1
        elif move == 'l':  # Move left
            if y > 0:
                if matrix[x][y-1] == matrix[x][y]:
                    matrix[x][y-1] *= 2
                    matrix[x][y] = 0
                elif matrix[x][y-1] == 0:
                    matrix[x][y-1] = matrix[x][y]
                    matrix[x][y] = 0
                y -= 1
        elif move == 'r':  # Move right
            if y < 3:
                if matrix[x][y+1] == matrix[x][y]:
                    matrix[x][y+1] *= 2
                    matrix[x][y] = 0
                elif matrix[x][y+1] == 0:
                    matrix[x][y+1] = matrix[x][y]
                    matrix[x][y] = 0
                y += 1
    return matrix

# Apply the movement sequence
final_matrix = move_number(matrix, position, sequence)

# Print the final matrix
print(final_matrix)