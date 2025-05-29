# Initial matrix
matrix = [
    [8, 4, 4],
    [0, 16, 0],
    [16, 0, 0]
]

# Initial position of the number to move
position = (1, 1)

# Sequence of moves
moves = 'llr'

# Function to perform the move
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
        if move == 'l':
            if y > 0:
                if matrix[x][y-1] == 0:
                    matrix[x][y-1] = matrix[x][y]
                    matrix[x][y] = 0
                elif matrix[x][y-1] == matrix[x][y]:
                    matrix[x][y-1] *= 2
                    matrix[x][y] = 0
                y -= 1
        elif move == 'r':
            if y < 2:
                if matrix[x][y+1] == 0:
                    matrix[x][y+1] = matrix[x][y]
                    matrix[x][y] = 0
                elif matrix[x][y+1] == matrix[x][y]:
                    matrix[x][y+1] *= 2
                    matrix[x][y] = 0
                y += 1
    return matrix

# Perform the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)