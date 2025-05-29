# Initial matrix
matrix = [
    [0, 0, 8],
    [32, 32, 32],
    [0, 0, 0]
]

# Initial position of the number
position = (1, 1)

# Move sequence
moves = 'lrl'

# Function to perform a move
def move_number(matrix, position, direction):
    x, y = position
    if direction == 'l':
        if y > 0:
            if matrix[x][y-1] == 0:
                matrix[x][y-1] = matrix[x][y]
                matrix[x][y] = 0
            elif matrix[x][y-1] == matrix[x][y]:
                matrix[x][y-1] *= 2
                matrix[x][y] = 0
            y -= 1
    elif direction == 'r':
        if y < 2:
            if matrix[x][y+1] == 0:
                matrix[x][y+1] = matrix[x][y]
                matrix[x][y] = 0
            elif matrix[x][y+1] == matrix[x][y]:
                matrix[x][y+1] *= 2
                matrix[x][y] = 0
            y += 1
    return (x, y)

# Perform the moves
for move in moves:
    position = move_number(matrix, position, move)

# Print the final matrix
print(matrix)