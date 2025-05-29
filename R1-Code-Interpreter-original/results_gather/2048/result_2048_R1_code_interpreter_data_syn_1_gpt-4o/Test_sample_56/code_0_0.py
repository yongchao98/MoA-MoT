# Initial matrix
matrix = [
    [0, 32, 8, 0],
    [0, 0, 0, 0],
    [8, 0, 0, 32],
    [16, 0, 0, 0]
]

# Initial position of the number to move
position = (2, 1)

# Sequence of moves
moves = 'dlul'

# Function to perform the move
def move_number(matrix, position, direction):
    x, y = position
    if direction == 'd':
        if x < 3 and matrix[x+1][y] == matrix[x][y]:
            matrix[x+1][y] *= 2
            matrix[x][y] = 0
            x += 1
        elif x < 3 and matrix[x+1][y] == 0:
            matrix[x+1][y] = matrix[x][y]
            matrix[x][y] = 0
            x += 1
    elif direction == 'l':
        if y > 0 and matrix[x][y-1] == matrix[x][y]:
            matrix[x][y-1] *= 2
            matrix[x][y] = 0
            y -= 1
        elif y > 0 and matrix[x][y-1] == 0:
            matrix[x][y-1] = matrix[x][y]
            matrix[x][y] = 0
            y -= 1
    elif direction == 'u':
        if x > 0 and matrix[x-1][y] == matrix[x][y]:
            matrix[x-1][y] *= 2
            matrix[x][y] = 0
            x -= 1
        elif x > 0 and matrix[x-1][y] == 0:
            matrix[x-1][y] = matrix[x][y]
            matrix[x][y] = 0
            x -= 1
    elif direction == 'r':
        if y < 3 and matrix[x][y+1] == matrix[x][y]:
            matrix[x][y+1] *= 2
            matrix[x][y] = 0
            y += 1
        elif y < 3 and matrix[x][y+1] == 0:
            matrix[x][y+1] = matrix[x][y]
            matrix[x][y] = 0
            y += 1
    return (x, y)

# Apply the sequence of moves
for move in moves:
    position = move_number(matrix, position, move)

# Print the final matrix
print(matrix)