# Initial matrix
matrix = [
    [0, 0, 4],
    [0, 2, 0],
    [0, 0, 16]
]

# Initial position of the number 2
position = (2, 1)

# Sequence of moves
moves = 'duu'

# Function to perform the move
def move_number(matrix, position, moves):
    for move in moves:
        x, y = position
        if move == 'u' and x > 0:
            # Move up
            if matrix[x-1][y] == matrix[x][y]:
                matrix[x-1][y] *= 2
                matrix[x][y] = 0
            elif matrix[x-1][y] == 0:
                matrix[x-1][y] = matrix[x][y]
                matrix[x][y] = 0
            position = (x-1, y)
        elif move == 'd' and x < 2:
            # Move down
            if matrix[x+1][y] == matrix[x][y]:
                matrix[x+1][y] *= 2
                matrix[x][y] = 0
            elif matrix[x+1][y] == 0:
                matrix[x+1][y] = matrix[x][y]
                matrix[x][y] = 0
            position = (x+1, y)
        elif move == 'l' and y > 0:
            # Move left
            if matrix[x][y-1] == matrix[x][y]:
                matrix[x][y-1] *= 2
                matrix[x][y] = 0
            elif matrix[x][y-1] == 0:
                matrix[x][y-1] = matrix[x][y]
                matrix[x][y] = 0
            position = (x, y-1)
        elif move == 'r' and y < 2:
            # Move right
            if matrix[x][y+1] == matrix[x][y]:
                matrix[x][y+1] *= 2
                matrix[x][y] = 0
            elif matrix[x][y+1] == 0:
                matrix[x][y+1] = matrix[x][y]
                matrix[x][y] = 0
            position = (x, y+1)
    return matrix

# Perform the moves
result_matrix = move_number(matrix, position, moves)

# Print the resulting matrix
print(result_matrix)