# Initial matrix
matrix = [
    [0, 0, 8, 0, 0],
    [0, 0, 32, 0, 32],
    [0, 32, 0, 0, 4],
    [0, 0, 0, 0, 4],
    [0, 0, 0, 32, 0]
]

# Initial position of the number to move
position = (4, 5)
# Sequence of moves
moves = 'lddud'

# Function to perform a move
def move_number(matrix, position, direction):
    x, y = position
    if direction == 'l':
        if y > 0 and (matrix[x][y-1] == 0 or matrix[x][y-1] == matrix[x][y]):
            matrix[x][y-1] += matrix[x][y]
            matrix[x][y] = 0
            y -= 1
    elif direction == 'r':
        if y < 4 and (matrix[x][y+1] == 0 or matrix[x][y+1] == matrix[x][y]):
            matrix[x][y+1] += matrix[x][y]
            matrix[x][y] = 0
            y += 1
    elif direction == 'u':
        if x > 0 and (matrix[x-1][y] == 0 or matrix[x-1][y] == matrix[x][y]):
            matrix[x-1][y] += matrix[x][y]
            matrix[x][y] = 0
            x -= 1
    elif direction == 'd':
        if x < 4 and (matrix[x+1][y] == 0 or matrix[x+1][y] == matrix[x][y]):
            matrix[x+1][y] += matrix[x][y]
            matrix[x][y] = 0
            x += 1
    return (x, y)

# Apply the sequence of moves
for move in moves:
    position = move_number(matrix, position, move)

# Print the final matrix
print(matrix)