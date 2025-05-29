# Initial matrix
matrix = [
    [0, 4, 0],
    [16, 16, 4],
    [0, 32, 0]
]

# Initial position of the number to move
x, y = 1, 1

# Move sequence
moves = 'ull'

# Function to perform a move
def move_number(matrix, x, y, direction):
    if direction == 'u' and x > 0:
        if matrix[x-1][y] == matrix[x][y]:
            matrix[x-1][y] *= 2
            matrix[x][y] = 0
        elif matrix[x-1][y] == 0:
            matrix[x-1][y] = matrix[x][y]
            matrix[x][y] = 0
        x -= 1
    elif direction == 'l' and y > 0:
        if matrix[x][y-1] == matrix[x][y]:
            matrix[x][y-1] *= 2
            matrix[x][y] = 0
        elif matrix[x][y-1] == 0:
            matrix[x][y-1] = matrix[x][y]
            matrix[x][y] = 0
        y -= 1
    return x, y

# Execute the moves
for move in moves:
    x, y = move_number(matrix, x, y, move)

# Print the final matrix
print(matrix)