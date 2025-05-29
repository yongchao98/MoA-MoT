# Initial matrix
matrix = [
    [4, 32, 0, 16, 0],
    [16, 0, 0, 16, 0],
    [16, 4, 16, 0, 0],
    [0, 0, 32, 16, 0],
    [0, 0, 0, 0, 0]
]

# Initial position of the number to move
x, y = 3, 3

# Sequence of moves
moves = 'rrdluu'

# Function to perform a move
def move(x, y, direction):
    if direction == 'r':
        if y < 4 and (matrix[x][y+1] == 0 or matrix[x][y+1] == matrix[x][y]):
            if matrix[x][y+1] == matrix[x][y]:
                matrix[x][y+1] *= 2
                matrix[x][y] = 0
            else:
                matrix[x][y+1] = matrix[x][y]
                matrix[x][y] = 0
            y += 1
    elif direction == 'l':
        if y > 0 and (matrix[x][y-1] == 0 or matrix[x][y-1] == matrix[x][y]):
            if matrix[x][y-1] == matrix[x][y]:
                matrix[x][y-1] *= 2
                matrix[x][y] = 0
            else:
                matrix[x][y-1] = matrix[x][y]
                matrix[x][y] = 0
            y -= 1
    elif direction == 'u':
        if x > 0 and (matrix[x-1][y] == 0 or matrix[x-1][y] == matrix[x][y]):
            if matrix[x-1][y] == matrix[x][y]:
                matrix[x-1][y] *= 2
                matrix[x][y] = 0
            else:
                matrix[x-1][y] = matrix[x][y]
                matrix[x][y] = 0
            x -= 1
    elif direction == 'd':
        if x < 4 and (matrix[x+1][y] == 0 or matrix[x+1][y] == matrix[x][y]):
            if matrix[x+1][y] == matrix[x][y]:
                matrix[x+1][y] *= 2
                matrix[x][y] = 0
            else:
                matrix[x+1][y] = matrix[x][y]
                matrix[x][y] = 0
            x += 1
    return x, y

# Perform the sequence of moves
for move_direction in moves:
    x, y = move(x, y, move_direction)

# Print the final matrix
print(matrix)