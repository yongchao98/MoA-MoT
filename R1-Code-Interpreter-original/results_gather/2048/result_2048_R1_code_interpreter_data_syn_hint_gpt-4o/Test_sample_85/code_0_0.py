# Initial matrix
matrix = [
    [0, 2, 0, 0, 4],
    [32, 16, 0, 2, 32],
    [0, 2, 0, 0, 0],
    [0, 2, 0, 0, 0],
    [8, 0, 2, 0, 0]
]

# Initial position
x, y = 2, 1

# Sequence of moves
moves = 'ddldlu'

# Function to perform the move
def move(x, y, direction):
    if direction == 'u' and x > 0:
        if matrix[x-1][y] == matrix[x][y]:
            matrix[x-1][y] *= 2
            matrix[x][y] = 0
        elif matrix[x-1][y] == 0:
            matrix[x-1][y] = matrix[x][y]
            matrix[x][y] = 0
        x -= 1
    elif direction == 'd' and x < 4:
        if matrix[x+1][y] == matrix[x][y]:
            matrix[x+1][y] *= 2
            matrix[x][y] = 0
        elif matrix[x+1][y] == 0:
            matrix[x+1][y] = matrix[x][y]
            matrix[x][y] = 0
        x += 1
    elif direction == 'l' and y > 0:
        if matrix[x][y-1] == matrix[x][y]:
            matrix[x][y-1] *= 2
            matrix[x][y] = 0
        elif matrix[x][y-1] == 0:
            matrix[x][y-1] = matrix[x][y]
            matrix[x][y] = 0
        y -= 1
    elif direction == 'r' and y < 4:
        if matrix[x][y+1] == matrix[x][y]:
            matrix[x][y+1] *= 2
            matrix[x][y] = 0
        elif matrix[x][y+1] == 0:
            matrix[x][y+1] = matrix[x][y]
            matrix[x][y] = 0
        y += 1
    return x, y

# Perform the moves
for move_direction in moves:
    x, y = move(x, y, move_direction)

# Print the final matrix
print(matrix)