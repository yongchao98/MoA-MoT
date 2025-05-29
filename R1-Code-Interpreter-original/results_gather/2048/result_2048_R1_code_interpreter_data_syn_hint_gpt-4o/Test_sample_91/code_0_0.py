# Initial matrix
matrix = [
    [8, 0, 0, 0, 32],
    [32, 0, 8, 0, 0],
    [0, 0, 0, 32, 2],
    [0, 32, 32, 0, 0],
    [0, 0, 2, 0, 0]
]

# Initial position of the number
position = (4, 2)  # 0-indexed, corresponds to (5, 3) in 1-indexed

# Sequence of moves
moves = 'rllddl'

# Function to perform a move
def move(position, direction):
    x, y = position
    if direction == 'r':
        if y < 4 and matrix[x][y + 1] == 0:
            matrix[x][y + 1] = matrix[x][y]
            matrix[x][y] = 0
            y += 1
    elif direction == 'l':
        if y > 0 and matrix[x][y - 1] == 0:
            matrix[x][y - 1] = matrix[x][y]
            matrix[x][y] = 0
            y -= 1
    elif direction == 'd':
        if x < 4 and matrix[x + 1][y] == 0:
            matrix[x + 1][y] = matrix[x][y]
            matrix[x][y] = 0
            x += 1
    elif direction == 'u':
        if x > 0 and matrix[x - 1][y] == 0:
            matrix[x - 1][y] = matrix[x][y]
            matrix[x][y] = 0
            x -= 1
    return (x, y)

# Perform the moves
for move_direction in moves:
    position = move(position, move_direction)

# Print the final matrix
print(matrix)