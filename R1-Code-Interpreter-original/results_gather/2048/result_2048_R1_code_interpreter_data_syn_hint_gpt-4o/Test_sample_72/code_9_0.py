# Initial matrix
matrix = [
    [0, 0, 4, 0, 0],
    [0, 0, 32, 0, 32],
    [0, 0, 0, 0, 0],
    [8, 0, 0, 0, 16],
    [0, 32, 0, 0, 4]
]

# Initial position of the number to move
position = (4, 3)

# Sequence of moves
moves = 'lluud'

# Function to perform the move
def move_number(matrix, position, moves):
    x, y = position
    for move in moves:
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

# Update the matrix based on the moves
final_matrix = move_number(matrix, position, moves)

# Print the final matrix
print(final_matrix)