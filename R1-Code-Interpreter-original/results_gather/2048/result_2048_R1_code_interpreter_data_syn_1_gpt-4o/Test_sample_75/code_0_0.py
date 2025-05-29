# Initial matrix
matrix = [
    [0, 16, 4, 16, 0],
    [0, 2, 0, 0, 32],
    [0, 0, 0, 0, 16],
    [4, 0, 0, 0, 8],
    [0, 0, 0, 2, 0]
]

# Initial position of the number
position = (3, 3)

# Move sequence
moves = 'urrld'

# Function to perform a move
def move(position, direction):
    x, y = position
    if direction == 'u':
        x -= 1
    elif direction == 'd':
        x += 1
    elif direction == 'l':
        y -= 1
    elif direction == 'r':
        y += 1
    return x, y

# Function to check if a move is valid
def is_valid(x, y):
    return 0 <= x < 5 and 0 <= y < 5

# Simulate the moves
for move_direction in moves:
    new_x, new_y = move(position, move_direction)
    if is_valid(new_x, new_y):
        # Check if we can combine
        if matrix[new_x][new_y] == matrix[position[0]][position[1]]:
            matrix[new_x][new_y] *= 2
            matrix[position[0]][position[1]] = 0
        elif matrix[new_x][new_y] == 0:
            matrix[new_x][new_y] = matrix[position[0]][position[1]]
            matrix[position[0]][position[1]] = 0
        position = (new_x, new_y)

# Print the final matrix
print(matrix)