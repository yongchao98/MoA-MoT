# Initial matrix
matrix = [
    [0, 0, 0, 16, 16],
    [0, 0, 0, 2, 2],
    [0, 0, 0, 0, 0],
    [0, 0, 16, 4, 0],
    [0, 0, 0, 0, 0]
]

# Starting position and value
x, y = 1, 3
value = matrix[x][y]

# Move sequence
moves = 'rlrdr'

# Function to perform a move
def move(x, y, direction):
    if direction == 'r':
        y += 1
    elif direction == 'l':
        y -= 1
    elif direction == 'd':
        x += 1
    elif direction == 'u':
        x -= 1
    return x, y

# Simulate the moves
for move_direction in moves:
    new_x, new_y = move(x, y, move_direction)
    if 0 <= new_x < 5 and 0 <= new_y < 5:
        if matrix[new_x][new_y] == value:
            value *= 2
            matrix[new_x][new_y] = 0
        x, y = new_x, new_y

# Place the final value in the matrix
matrix[x][y] = value

# Print the final matrix
print(matrix)