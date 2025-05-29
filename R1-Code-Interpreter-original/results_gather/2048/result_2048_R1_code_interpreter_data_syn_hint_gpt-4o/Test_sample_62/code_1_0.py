# Initial matrix
matrix = [
    [8, 0, 0, 32, 0],
    [4, 0, 0, 16, 0],
    [0, 0, 0, 0, 0],
    [2, 4, 0, 2, 4],
    [4, 0, 0, 0, 0]
]

# Initial position of the number to move
position = (2, 4)

# Sequence of moves
moves = 'dldur'

# Function to move the position
def move(position, direction):
    x, y = position
    if direction == 'u':
        return (max(x - 1, 0), y)
    elif direction == 'd':
        return (min(x + 1, 4), y)
    elif direction == 'l':
        return (x, max(y - 1, 0))
    elif direction == 'r':
        return (x, min(y + 1, 4))
    return position

# Function to perform a move and merge if possible
def perform_move(matrix, position, direction):
    x, y = position
    new_x, new_y = move(position, direction)
    if matrix[new_x][new_y] == matrix[x][y]:
        matrix[new_x][new_y] *= 2
        matrix[x][y] = 0
    elif matrix[new_x][new_y] == 0:
        matrix[new_x][new_y] = matrix[x][y]
        matrix[x][y] = 0
    return (new_x, new_y)

# Simulate the moves
for move_direction in moves:
    position = perform_move(matrix, position, move_direction)

# Print the final matrix
print(matrix)