# Initial matrix
matrix = [
    [8, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 8, 0],
    [0, 0, 0, 0, 0],
    [0, 32, 0, 16, 32]
]

# Initial position of the number to move
x, y = 4, 3

# Sequence of moves
moves = 'uulll'

# Function to move the number
def move_number(matrix, x, y, moves):
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
            if x < len(matrix) - 1:
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
            if y < len(matrix[0]) - 1:
                if matrix[x][y+1] == matrix[x][y]:
                    matrix[x][y+1] *= 2
                    matrix[x][y] = 0
                elif matrix[x][y+1] == 0:
                    matrix[x][y+1] = matrix[x][y]
                    matrix[x][y] = 0
                y += 1
    return matrix

# Move the number according to the sequence
result_matrix = move_number(matrix, x, y, moves)

# Print the resulting matrix
print(result_matrix)