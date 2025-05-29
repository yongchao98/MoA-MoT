# Initial matrix
matrix = [
    [0, 0, 0, 8, 0],
    [0, 16, 0, 16, 0],
    [0, 0, 0, 0, 8],
    [0, 2, 2, 0, 16],
    [2, 0, 0, 2, 0]
]

# Initial position
x, y = 3, 3

# Movement sequence
sequence = 'rrldd'

# Function to move and combine numbers
def move_and_combine(matrix, x, y, direction):
    if direction == 'r':
        if y < 4:
            if matrix[x][y] == matrix[x][y + 1]:
                matrix[x][y + 1] *= 2
                matrix[x][y] = 0
            elif matrix[x][y + 1] == 0:
                matrix[x][y + 1] = matrix[x][y]
                matrix[x][y] = 0
            y += 1
    elif direction == 'l':
        if y > 0:
            if matrix[x][y] == matrix[x][y - 1]:
                matrix[x][y - 1] *= 2
                matrix[x][y] = 0
            elif matrix[x][y - 1] == 0:
                matrix[x][y - 1] = matrix[x][y]
                matrix[x][y] = 0
            y -= 1
    elif direction == 'd':
        if x < 4:
            if matrix[x][y] == matrix[x + 1][y]:
                matrix[x + 1][y] *= 2
                matrix[x][y] = 0
            elif matrix[x + 1][y] == 0:
                matrix[x + 1][y] = matrix[x][y]
                matrix[x][y] = 0
            x += 1
    elif direction == 'u':
        if x > 0:
            if matrix[x][y] == matrix[x - 1][y]:
                matrix[x - 1][y] *= 2
                matrix[x][y] = 0
            elif matrix[x - 1][y] == 0:
                matrix[x - 1][y] = matrix[x][y]
                matrix[x][y] = 0
            x -= 1
    return x, y

# Execute the movement sequence
for move in sequence:
    x, y = move_and_combine(matrix, x, y, move)

# Print the final matrix
print(matrix)