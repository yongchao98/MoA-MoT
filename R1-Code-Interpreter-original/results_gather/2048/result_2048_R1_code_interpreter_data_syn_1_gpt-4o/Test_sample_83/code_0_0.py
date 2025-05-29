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
moves = 'lullll'

# Function to move the number in the matrix
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'u':
            if x > 0:
                if matrix[x-1][y] == matrix[x][y]:
                    matrix[x-1][y] *= 2
                    matrix[x][y] = 0
                elif matrix[x-1][y] == 0:
                    matrix[x-1][y] = matrix[x][y]
                    matrix[x][y] = 0
                x -= 1
        elif move == 'd':
            if x < 4:
                if matrix[x+1][y] == matrix[x][y]:
                    matrix[x+1][y] *= 2
                    matrix[x][y] = 0
                elif matrix[x+1][y] == 0:
                    matrix[x+1][y] = matrix[x][y]
                    matrix[x][y] = 0
                x += 1
        elif move == 'l':
            if y > 0:
                if matrix[x][y-1] == matrix[x][y]:
                    matrix[x][y-1] *= 2
                    matrix[x][y] = 0
                elif matrix[x][y-1] == 0:
                    matrix[x][y-1] = matrix[x][y]
                    matrix[x][y] = 0
                y -= 1
        elif move == 'r':
            if y < 4:
                if matrix[x][y+1] == matrix[x][y]:
                    matrix[x][y+1] *= 2
                    matrix[x][y] = 0
                elif matrix[x][y+1] == 0:
                    matrix[x][y+1] = matrix[x][y]
                    matrix[x][y] = 0
                y += 1
    return matrix

# Apply the moves
final_matrix = move_number(matrix, position, moves)

# Print the final matrix
print(final_matrix)